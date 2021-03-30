classdef Electrode < PhysicalModel
    
    properties
        
        %% submodels:
        % electrode active component
        ElectrodeActiveComponent %  (class ActiveElectroChemicalComponent)
        % current collector
        CurrentCollector % (class ElectronicComponent)

        couplingTerms;
        
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
            
            eac = model.ElectrodeActiveComponent;
            cc = model.CurrentCollector;            

            eac_phi = state.ElectrodeActiveComponent.phi;
            cc_phi = state.CurrentCollector.phi;

            eac_sigmaeff = ne.EffectiveElectronicConductivity;
            cc_sigmaeff = cc.EffectiveElectronicConductivity;
    
            %% We setup the current transfers between CurrentCollector and ElectrodeActiveComponent
            
            eac_j_bcsource  = eac_phi*0.0; %NB hack to initialize zero ad
            cc_j_bcsource = cc_phi*0.0; %NB hack to initialize zero ad

            coupterm = getCoupTerm(model, 'CurrentCollector-ElectrodeActiveComponent');
            face_cc = coupterm.couplingfaces(:, 1);
            face_eac = coupterm.couplingfaces(:, 2);
            [teac, bccell_eac] = eac.operators.harmFaceBC(eac_sigmaeff, face_eac);
            [tcc, bccell_cc] = cc.operators.harmFaceBC(cc_sigmaeff, face_cc);

            bcphi_eac = eac_phi(bccell_eac);
            bcphi_cc = cc_phi(bccell_cc);

            trans = 1./(1./teac + 1./tcc);
            crosscurrent = trans.*(bcphi_cc - bcphi_eac);
            eac_j_bcsource(bccell_eac) = crosscurrent;
            cc_j_bcsource(bccell_cc) = -crosscurrent;


            state.ElectrodeActiveComponent.jBcSource = state.ElectrodeActiveComponent.jBcSource + eac_j_bcsource;
            state.CurrentCollector.jBcSource = state.CurrentCollector.jBcSource + cc_j_bcsource;
            
        end
        
        function state = updateT(model, state)
            names = {'ElectrodeActiveComponent', 'CurrentCollector'};
            for ind = 1 : numel(names)
                ame = names{ind};
                nc = model.(name).G.cells.num;
                state.(name).T = state.T(1)*ones(nc, 1);
            end
        end        
        
    end    
end

