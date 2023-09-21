classdef ProtonicMembraneAnode < ProtonicMembraneElectrode
    
    properties
        
        % coefficient in Buttler-Volmer
        beta
        
        % charge-transfer current densigy
        i0

        % Limiting current densities
        ila % Anode
        ilc % Cathode

        Rct
        ila0
        n
        SU
        pH2O_in
        
        
    end
    
    methods
        
        function model = ProtonicMembraneAnode(paramobj)

            model = model@ProtonicMembraneElectrode(paramobj);

            fdnames = {'beta', ...
                       'ila0', ...
                       'ilc' , ...
                       'Rct' , ...
                       'n'   , ...
                       'SU'  , ...
                       'pH2O_in'};
            model = dispatchParams(model, paramobj, fdnames);
            
            model.constants = PhysicalConstants();

            con = model.constants;

            f = con.F/(con.R*model.T);
            
            % Compute i0
            
            i0 = 1./(model.Rct.*f.*model.n); % charge-transfer current density

            % Compute ila
            
            pH2O = model.pH2O_in*(1 - model.SU);
            ila = model.ila0*pH2O;

            % Assign the values
            
            model.i0   = i0;
            model.ila  = ila;
            
        end
        
        function model = registerVarAndPropfuncNames(model)
        
            model = registerVarAndPropfuncNames@ProtonicMembraneElectrode(model);

            varnames = {'alpha'};

            model = model.registerVarNames(varnames);
            
            fn = @ProtonicMembraneAnode.updateJHp;
            inputnames = {'eta', 'alpha'};
            model = model.registerPropFunction({'jHp', fn, inputnames});
            
            % fn = @ProtonicMembraneAnode.updateEta2;
            % inputnames = {'j', 'alpha'};
            % model = model.registerPropFunction({'eta', fn, inputnames});

            % fn = @ProtonicMembraneAnode.updatePhi2;
            % inputnames = {'eta', 'pi', 'Eocp'};
            % model = model.registerPropFunction({'phi', fn, inputnames});


        end
        
        
        function state = updateJHp(model, state)
            
            con  = model.constants;

            beta = model.beta;
            ila  = model.ila;
            ilc  = model.ilc;
            i0   = model.i0;
            
            eta   = state.eta;
            alpha = state.alpha;
            
            feta = con.F*model.n/(con.R*model.T).*eta;
            
            jHp = -i0*(exp(-beta*feta) - exp((1 - beta)*feta))./(1 + (i0/ilc)*exp(-beta*feta) - (i0/ila)*exp(-(1 - beta)*feta));

            % jHp = i0*(exp(feta/2) - exp(-feta/2))/2;
            
            % R = 5;
            % jHp = 1/R*eta;
            
            state.jHp = jHp;
            
        end

        function state = updateEta2(model, state)

            R = 0.05;
            
            j = state.j;
            eta = R*j; 
            
            state.eta = eta;
            
        end

        function state = updatePhi2(model, state)

            state.phi = state.pi - state.eta - state.Eocp;
            
        end

        
    end
    
end

