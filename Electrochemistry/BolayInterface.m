classdef BolayInterface < Interface 

    properties

        SEImolarVolume
        SEIionicConductivity
        SEIelectronicDiffusionCoefficient
        SEIintersticialConcentration
        SEIstoichiometricCoefficient
        SEIlengthInitial
        
        SEIlengthRef
        SEIvoltageDropRef

    end

    methods

        function model = BolayInterface(inputparams)

            model = model@Interface(inputparams);

            fdnames = {'SEImolarVolume'                   , ...
                       'SEIionicConductivity'             , ...
                       'SEIelectronicDiffusionCoefficient', ...        
                       'SEIintersticialConcentration'     , ...
                       'SEIstoichiometricCoefficient'       , ...
                       'SEIlengthInitial'};
            
            model = dispatchParams(model, inputparams, fdnames);

            L0 = model.SEIlengthInitial;
            
            model.SEIlengthRef = L0;
            
        end

        function model = registerVarAndPropfuncNames(model)

            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@Interface(model);

            varnames = {};
            % Length of SEI layer [m]
            varnames{end + 1} = 'SEIlength';
            % normalized length of SEI layer []
            varnames{end + 1} = 'normalizedSEIlength';
            % potential drop at SEI [V]
            varnames{end + 1} = 'SEIvoltageDrop';
            % potential drop at SEI []
            varnames{end + 1} = 'normalizedSEIvoltageDrop';
            % SEI flux [mol/m^2/s]
            varnames{end + 1} = 'SEIflux';
            % SEI mass conservation
            varnames{end + 1} = 'SEImassCons';
            % potential in electrolyte
            varnames{end + 1} = 'SEIvoltageDropEquation';
            % concentration of Lithium in the SEI (per total volume) in mol/m^3
            varnames{end + 1} = 'SEIconcentration';
            
            model = model.registerVarNames(varnames);

            fn = @BolayInterface.updateSEIlength;
            inputnames = {'normalizedSEIlength'};
            model = model.registerPropFunction({'SEIlength', fn, inputnames});

            fn = @BolayInterface.updateSEIvoltageDrop;
            inputnames = {'normalizedSEIvoltageDrop'};
            model = model.registerPropFunction({'SEIvoltageDrop', fn, inputnames});
            
            fn = @BolayInterface.updateSEIflux;
            inputnames = {'SEIlength', 'SEIvoltageDrop', 'phiElectrode', 'phiElectrolyte', 'T'};
            model = model.registerPropFunction({'SEIflux', fn, inputnames});

            fn = @BolayInterface.updateSEImassCons;
            fn = {fn, @(propfunction) PropFunction.accumFuncCallSetupFn(propfunction)};
            inputnames = {'SEIflux', 'SEIlength'};
            model = model.registerPropFunction({'SEImassCons', fn, inputnames});

            fn = @BolayInterface.updateSEIvoltageDropEquation;
            inputnames = {'intercalationFlux', 'SEIvoltageDrop'};
            model = model.registerPropFunction({'SEIvoltageDropEquation', fn, inputnames});

            fn = @BolayInterface.updateEta;
            inputnames = {'phiElectrolyte', 'phiElectrode', 'OCP', 'SEIvoltageDrop'};
            model = model.registerPropFunction({'eta', fn, inputnames});

            fn = @BolayInterface.updateSEIconcentration;
            model = model.registerPropFunction({'SEIconcentration', fn, {'SEIlength'}});

            model = model.setAsExtraVarName('SEIconcentration');
            
        end


        function state = updateSEIconcentration(model, state)
            
            vsa     = model.volumetricSurfaceArea;
	    scoef   = model.SEIstoichiometricCoefficient;
	    seimvol = model.SEImolarVolume;

            l = state.SEIlength;
            
            state.SEIconcentration = (scoef/seimvol)*l.*vsa;

        end
        
        function state = updateSEIflux(model, state)

            De  = model.SEIelectronicDiffusionCoefficient;
            ce0 = model.SEIintersticialConcentration;

            R = model.constants.R;
            F = model.constants.F;            

            T              = state.T;
            phiElectrode   = state.phiElectrode;
            phiElectrolyte = state.phiElectrolyte;
            Usei           = state.SEIvoltageDrop;
            L              = state.SEIlength;

            eta = phiElectrode - phiElectrolyte - Usei;

            state.SEIflux = De*ce0./L.*exp(-(F./(R*T)).*eta).*(1 - (F./(2*R*T)).*Usei);
           
        end

        function state = updateSEIlength(model, state)

            state.SEIlength = model.SEIlengthRef*state.normalizedSEIlength;
            
        end

        function state = updateSEIvoltageDrop(model, state)

            state.SEIvoltageDrop = model.SEIvoltageDropRef*state.normalizedSEIvoltageDrop;
            
        end
        
        function state = updateSEImassCons(model, state, state0, dt)

            s = model.SEIstoichiometricCoefficient;
            V = model.SEImolarVolume;

            L0 = state0.SEIlength;
            
            L  = state.SEIlength;
            N  = state.SEIflux;

            state.SEImassCons = s/V.*(L - L0)./dt - N;
            
        end

        function newstate = addVariablesAfterConvergence(model, newstate, state)

            newstate = addVariablesAfterConvergence@Interface(model, newstate, state);
            newstate.SEIlength = state.SEIlength;
            
        end
            
        function state = updateSEIvoltageDropEquation(model, state)

            k = model.SEIionicConductivity;
            F = model.constants.F;

            U = state.SEIvoltageDrop;
            L = state.SEIlength;
            R = state.intercalationFlux;

            state.SEIvoltageDropEquation = U - F*R.*L/k;
            
        end

        function state = updateEta(model, state)

            phiElyte       = state.phiElectrolyte;
            phiElde        = state.phiElectrode;
            OCP            = state.OCP;
            SEIvoltageDrop = state.SEIvoltageDrop;

            state.eta = (phiElde - phiElyte - OCP - SEIvoltageDrop);
            
        end

    
    end
end

