classdef Electrode < PhysicalModel
    
    properties
        
        %% submodels:
        % electrode active component
        ElectrodeActiveComponent %  (class ElectrodeActiveComponent)
        % current collector
        CurrentCollector % (class ElectronicComponent)

        couplingTerm;
        
    end

    methods
        
        function model = Electrode(paramobj)
        % paramobj is instance of ElectrodeInputParams
            model = model@PhysicalModel([]);
            
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks', false);
            
            fdnames = {'G', ...
                       'couplingTerm'};
            model = dispatchParams(model, paramobj, fdnames);
            
            % Assign the two components
            model.ElectrodeActiveComponent = model.setupElectrodeActiveComponent(paramobj.eac);
            model.CurrentCollector = model.setupCurrentCollector(paramobj.cc);
           
        end
        
        function eac = setupElectrodeActiveComponent(model, paramobj)
        % paramobj is instanceo of ElectrodeActiveComponentInputParams
        % standard instantiation (ActiveMaterial is specified in ElectrodeActiveComponent instantiation)
            eac = ElectrodeActiveComponent(paramobj);
        end
        
        function cc = setupCurrentCollector(model, paramobj)
        % standard instantiation 
            cc = CurrentCollector(paramobj);
        end
        
        function state = setupCoupling(model, state)
        % setup coupling terms between the current collector and the electrode active component            
            
            elde  = model;
            eac   = 'ElectrodeActiveComponent';
            cc    = 'CurrentCollector';

            eac_phi = state.(eac).phi;
            cc_phi = state.(cc).phi;

            eac_sigmaeff = elde.(eac).EffectiveElectronicConductivity;
            cc_sigmaeff = elde.(cc).EffectiveElectronicConductivity;
    
            %% We setup the current transfers between CurrentCollector and ElectrodeActiveComponent
            
            eac_jCoupling  = eac_phi*0.0; %NB hack to initialize zero ad
            cc_jCoupling = cc_phi*0.0; %NB hack to initialize zero ad

            coupterm = model.couplingTerm;
            face_cc = coupterm.couplingfaces(:, 1);
            face_eac = coupterm.couplingfaces(:, 2);
            [teac, bccell_eac] = elde.(eac).operators.harmFaceBC(eac_sigmaeff, face_eac);
            [tcc, bccell_cc] = elde.(cc).operators.harmFaceBC(cc_sigmaeff, face_cc);

            bcphi_eac = eac_phi(bccell_eac);
            bcphi_cc = cc_phi(bccell_cc);

            trans = 1./(1./teac + 1./tcc);
            crosscurrent = trans.*(bcphi_cc - bcphi_eac);
            eac_jCoupling(bccell_eac) = crosscurrent;
            cc_jCoupling(bccell_cc) = -crosscurrent;

            % We set here volumetric current source to zero for current collector (could have been done at a more logical place but
            % let us do it here, for simplicity)
            state.(cc).eSource = zeros(elde.(cc).G.cells.num, 1);
            
            state.(eac).jCoupling = eac_jCoupling;
            state.(cc).jCoupling  = cc_jCoupling;

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

