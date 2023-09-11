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
                       'pH20_in'};
            model = dispatchParams(model, paramobj, fdnames);
            
            model.constants = PhysicalConstants();

            con = model.constants;

            f = con.F/(con.R*model.T);
            
            % Compute i0
            
            i0 = 1./(model.Rct.*f.*model.n); % charge-transfer current density

            % Compute ila
            
            pH2O = model.pH2O_in*(1 - model.SU);
            ila = model.ila0*pH2O;
            
            
        end
        
        function model = registerVarAndPropfuncNames(model)
        
            model = registerVarAndPropfuncNames@ProtonicMembraneElectrode(model);
            
            fn = @ProtonicMembraneAnode.updateJHp;
            inputnames = {'eta'};
            model = model.registerPropFunction({'jHp', fn, inputnames});

        end
        
        
        function state = updateJHp(model, state)
            
            con  = model.constants;

            beta = model.beta;
            ila  = model.ila;
            ilc  = model.ilc;
            i0   = model.i0;
            
            eta = state.eta;
            
            feta = con.F*model.n/(con.R*model.T).*eta;
            
            jHo = i0*(exp(-beta*feta) - exp((1 - beta)*feta))/(1 + (i0/ilc)*exp(-beta*feta) - (i0/ila)*exp(-(1 - beta)*feta));
            
            state.jHp = jHp;
            
        end
        
    end
    
end
