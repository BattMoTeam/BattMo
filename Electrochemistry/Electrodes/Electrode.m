classdef Electrode < PhysicalModel
    
    properties
        
        %% submodels:
        % electrode active component
        ElectrodeActiveComponent %  (class ActiveElectroChemicalComponent)
        % current collector
        CurrentCollector % (class ElectronicComponent)

        couplingTerms;
        couplingNames;
        
    end

    methods
        
        function model = Electrode(paramobj)
            
            model = model@PhysicalModel([]);
            
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks', true);
            
            fdnames = {'G', ...
                       'couplingTerms'};
            model = dispatchParams(model, paramobj, fdnames);
            
            % Assign the two components
            model.ElectrodeActiveComponent = ActiveElectroChemicalComponent(paramobj.eac);
            model.CurrentCollector = CurrentCollector(paramobj.cc);
            
            % setup couplingNames
            model.couplingNames = cellfun(@(x) x.name, model.couplingTerms, 'uniformoutput', false);
           
        end
        
        
        function state = setupCoupling(model, state)
        % initiate jBcSource

            elde  = model;
            eac   = 'ElectrodeActiveComponent';
            cc    = 'CurrentCollector';
            
            eac_jBcSource = zeros(elde.(eac).G.cells.num, 1);
            cc_jBcSource = zeros(elde.(cc).G.cells.num, 1);
            cc_eSource = zeros(elde.(cc).G.cells.num, 1);
            
            phi = state.(cc).phi;
            if isa(phi, 'ADI')
                adsample = getSampleAD(phi);
                adbackend = model.AutoDiffBackend;
                eac_jBcSource = adbackend.convertToAD(eac_jBcSource, adsample);
                cc_jBcSource = adbackend.convertToAD(cc_jBcSource, adsample);
            end
            
            state.(eac).jBcSource = eac_jBcSource;
            state.(cc).jBcSource = cc_jBcSource;
            state.(cc).eSource = cc_eSource; % always zero, no volumetric electron source
            
            state = model.setupCurrentCollectorCoupling(state);
            state = model.setupBoundaryCoupling(state);
        end

        function state = setupBoundaryCoupling(model, state);
        % We impose the boundary condition at chosen boundary cells of the current collector
        % shortcuts:
        % cc  : CurrentCollector            
            
            coupnames = model.couplingNames;
            coupterm = getCoupTerm(model.couplingTerms, 'bc-CurrentCollector', coupnames);
            
            elde = model;
            cc = 'CurrentCollector';
            
            phi = state.(cc).phi;
            
            j_bcsource = phi*0.0; %NB hack to initialize zero ad
            
            sigmaeff = elde.(cc).EffectiveElectronicConductivity;
            faces = coupterm.couplingfaces;
            bcval = zeros(numel(faces), 1);
            [t, cells] = elde.(cc).operators.harmFaceBC(sigmaeff, faces);
            j_bcsource(cells) = j_bcsource(cells) + t.*(bcval - phi(cells));
            
            state.(cc).jBcSource = state.(cc).jBcSource + j_bcsource;
            
        end
        
        
        function state  = setupCurrentCollectorCoupling(model, state)
        % setup electrical coupling between current collector and the electrode active component
        % shortcuts:
        % eac : ElectrodeActiveComponent
        % cc  : CurrentCollector
            
            elde  = model;
            eac   = 'ElectrodeActiveComponent';
            cc    = 'CurrentCollector';

            eac_phi = state.(eac).phi;
            cc_phi = state.(cc).phi;

            eac_sigmaeff = elde.(eac).EffectiveElectronicConductivity;
            cc_sigmaeff = elde.(cc).EffectiveElectronicConductivity;
    
            %% We setup the current transfers between CurrentCollector and ElectrodeActiveComponent
            
            eac_j_bcsource  = eac_phi*0.0; %NB hack to initialize zero ad
            cc_j_bcsource = cc_phi*0.0; %NB hack to initialize zero ad

            coupnames = model.couplingNames;
            coupterm = getCoupTerm(model.couplingTerms, 'CurrentCollector-ElectrodeActiveComponent', coupnames);
            face_cc = coupterm.couplingfaces(:, 1);
            face_eac = coupterm.couplingfaces(:, 2);
            [teac, bccell_eac] = elde.(eac).operators.harmFaceBC(eac_sigmaeff, face_eac);
            [tcc, bccell_cc] = elde.(cc).operators.harmFaceBC(cc_sigmaeff, face_cc);

            bcphi_eac = eac_phi(bccell_eac);
            bcphi_cc = cc_phi(bccell_cc);

            trans = 1./(1./teac + 1./tcc);
            crosscurrent = trans.*(bcphi_cc - bcphi_eac);
            eac_j_bcsource(bccell_eac) = crosscurrent;
            cc_j_bcsource(bccell_cc) = -crosscurrent;

            state.(eac).jBcSource = state.(eac).jBcSource + eac_j_bcsource;
            state.(cc).jBcSource = state.(cc).jBcSource + cc_j_bcsource;
            
        end
        
        function state = updateT(model, state)
            names = {'ElectrodeActiveComponent', 'CurrentCollector'};
            for ind = 1 : numel(names)
                name = names{ind};
                nc = model.(name).G.cells.num;
                state.(name).T = state.T(1)*ones(nc, 1);
            end
        end        
        
    end    
end

