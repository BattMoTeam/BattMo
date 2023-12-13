classdef MagnesiumElectrolyte < SeaWaterElectrolyte

    properties
        
        % K is given by
        % K(1) : 'H2O         = H+ + OH-'
        % K(2) : 'Mg+2 + Cl-  = MgCl+'
        % K(3) : 'Mg+2 + OH-  = MgOH+'
        % K(4) : 'Mg(OH)2     = Mg+2 + 2OH-'
    
    end
    
    methods

        function model = MagnesiumElectrolyte(inputparams)
            model = model@SeaWaterElectrolyte(inputparams);
        end

        function model = registerVarAndPropfuncNames(model)
            
            model = registerVarAndPropfuncNames@SeaWaterElectrolyte(model);

            spdict     = model.spdict;
            nsp        = model.nsp;
            nqp        = model.nqp;
            indmainsp  = model.indmainsp;
            indsolidsp = model.indsolidsp;
            
            %% update saturation concentration and reaction coefficient
            fn = @() MagnesiumElectrolyte.updatePrecipitationCoef;
            inputnames = {VarName({}, 'cs', nsp, [spdict('OH-'), spdict('Mg+2')])};
            model = model.registerPropFunction({{'k'}, fn, inputnames});
            model = model.registerPropFunction({{'cSat'}, fn, inputnames});            
            
            
            %% update the concentrations in the model using chemical constants
            fn = @() MagnesiumElectrolyte.updateConcentrations;
            indinput = [indmainsp, indsolidsp];
            inputnames = {VarName({}, 'cs', nsp, indinput)};
            indsec = setdiff((1 : nsp), indinput);
            model = model.registerPropFunction({VarName({}, 'cs', nsp, indsec), fn, inputnames});

        end

        function model = setupMainIonIndex(model)
            model.mainIonIndex = model.spdict('Mg+2');
        end
        
        function state = updatePrecipitationCoef(model, state)

            d    = model.spdict;
            K    = model.K(4);
            
            cs = state.cs;
            
            cOH_m  = cs{d('OH-')};
            cMg_p2 = cs{d('Mg+2')};
            
            cSat = K./(cOH_m.^2);
            
            % TODO : Hard-coded constants. Should be moved outside the function.
            D     = 1e-9;
            delta = 100e-6;
            
            state.k = D/delta*(cMg_p2 - cSat);
            state.cSat = cSat;
            
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
