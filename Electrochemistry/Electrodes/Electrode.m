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
            
            fdnames = {'G', ...
                       'couplingTerms'};
            model = dispatchParams(model, paramobj, fdnames);
            
            % Assign the two components
            model.ElectrodeActiveComponent = ActiveElectroChemicalComponent(paramobj.eac);
            model.CurrentCollector = ElectronicComponent(paramobj.cc);
            
            % setup couplingNames
            model.couplingNames = cellfun(@(x) x.name, model.couplingTerms, 'uniformoutput', false);
           
        end
        
        
        function state = setupCoupling(model, state)
            state = model.setupCurrentCollectorCoupling(state);
            state = model.setupBoundaryCoupling(state);
        end

        function state = setupBoundaryCoupling(model, state);
        % We impose the boundary condition at chosen boundary cells of the current collector
        % shortcuts:
        % cc  : CurrentCollector            
            
            coupterm = model.getCoupTerm('bc-NegativeCurrentCollector');
            
            cc = model.CurrentCollector;            
            phi = state.CurrentCollector.phi;
            
            j_bcsource = phi*0.0; %NB hack to initialize zero ad
            
            sigmaeff = cc.EffectiveElectronicConductivity;
            faces = coupterm.couplingfaces;
            bcval = zeros(numel(faces), 1);
            [t, cells] = cc.operators.harmFaceBC(sigmaeff, faces);
            j_bcsource(cells) = j_bcsource(cells) + t.*(bcval - phi(cells));
            
            state.CurrentCollector.jBcSource = state.CurrentCollector.jBcSource + cc_j_bcsource;
            
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

            state.(elde).(eac).jBcSource = state.(elde).(eac).jBcSource + eac_j_bcsource;
            state.(elde).(cc).jBcSource = state.(elde).(ccc).jBcSource + cc_j_bcsource;
            
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

