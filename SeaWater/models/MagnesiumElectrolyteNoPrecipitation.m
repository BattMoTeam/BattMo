classdef MagnesiumElectrolyteNoPrecipitation < SeaWaterElectrolyteNoPrecipitation

    properties
        
        % K is given by
        % K(1) : 'H2O         = H+ + OH-'
        % K(2) : 'Mg+2 + Cl-  = MgCl+'
        % K(3) : 'Mg+2 + OH-  = MgOH+'
        % K(4) : 'Mg(OH)2     = Mg+2 + 2OH-'
    
    end
    
    methods

        function model = MagnesiumElectrolyteNoPrecipitation(inputparams)
            model = model@SeaWaterElectrolyteNoPrecipitation(inputparams);
        end

        function model = registerVarAndPropfuncNames(model)
            
            model = registerVarAndPropfuncNames@SeaWaterElectrolyteNoPrecipitation(model);

            spdict     = model.spdict;
            nsp        = model.nsp;
            nqp        = model.nqp;
            indmainsp  = model.indmainsp;
            indsolidsp = model.indsolidsp;
            
            %% update the concentrations in the model using chemical constants
            fn = @() MagnesiumElectrolyteNoPrecipitation.updateConcentrations;
            indinput = [indmainsp, indsolidsp];
            inputnames = {VarName({}, 'cs', nsp, indinput)};
            indsec = setdiff((1 : nsp), indinput);
            model = model.registerPropFunction({VarName({}, 'cs', nsp, indsec), fn, inputnames});

        end

        function model = setupMainIonIndex(model)
            model.mainIonIndex = model.spdict('Mg+2');
        end
        
        function state = updateConcentrations(model, state)

            d = model.spdict; % dictionary with index of the species
            K = model.K;
            
            cs = state.cs;
            
            cH_p   = cs{d('H+')};
            cCl_m  = cs{d('Cl-')};
            cMg_p2 = cs{d('Mg+2')};
            
            cOH_m = K(1)./cH_p;
            
            cs{d('OH-')}   = cOH_m;
            cs{d('MgCl+')} = K(2).*cCl_m.*cMg_p2;
            cs{d('MgOH+')} = K(3).*cOH_m.*cMg_p2;

            state.cs = cs;
            
        end
        
    end
end
