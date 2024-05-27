classdef BolayInterface < Interface 

    properties

        SEImolarVolume
        SEIionicConductivity
        SEIelectronicDiffusionCoefficient
        SEIintersticialConcentration
        SEIstochiometricCoeffcient

    end

    methods

        function model = BolayInterface(inputparams)

            model = model@Interface(inputparams);

            fdnames = {'SEImolarVolume'                   , ...
                       'SEIionicConductivity'             , ...
                       'SEIelectronicDiffusionCoefficient', ...        
                       'SEIintersticialConcentration'     , ...
                       'SEIstochiometricCoeffcient'};
            
            model = dispatchParams(model, inputparams, fdnames);

        end

        function model = registerVarAndPropfuncNames(model)

            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@Interface(model);

            varnames = {};
            % Length of SEI layer [m]
            varnames{end + 1} = 'SEIlength';
            % potential drop at SEI [V]
            varnames{end + 1} = 'SEIvoltageDrop';
            % SEI flux [mol/m^2/s]
            varnames{end + 1} = 'SEIflux';
            % SEI mass conservation
            varnames{end + 1} = 'SEImassCons';
            % potential in electrolyte
            varnames{end + 1} = 'SEIvoltageDropEquation';

            model = model.registerVarNames(varnames);

            fn = @BolayInterface.updateSEIflux;
            inputnames = {'SEIlength', 'SEIvoltageDrop', 'eta'};
            model = model.registerPropFunction({'SEIflux', fn, inputnames});

            fn = @BolayInterface.updateSEImassCons;
            fn = {fn, @(propfunction) PropFunction.accumFuncCallSetupFn(propfunction)};
            inputnames = {'SEIflux', 'SEIlength'};
            model = model.registerPropFunction({'SEImassCons', fn, inputnames});

            fn = @BolayInterface.updateSEIvoltageDropEquation;
            inputnames = {'R', 'SEIvoltageDrop'};
            model = model.registerPropFunction({'SEIvoltageDropEquation', fn, inputnames});

            fn = @BolayInterface.updateEta;
            inputnames = {'phiElectrolyte', 'phiElectrode', 'OCP', 'SEIvoltageDrop'};
            model = model.registerPropFunction({'eta', fn, inputnames});

        end

        function state = updateSEIflux(model, state)

            De  = model.SEIelectronicDiffusionCoefficient;
            ce0 = model.SEIintersticialConcentration;

            R = model.constants.R;
            F = model.constants.F;            

            T   = state.T;
            eta = state.eta;
            U   = state.SEIvoltageDrop;
            L   = state.SEIlength;

            state.SEIflux = De*ce0./L.*exp(-(F./(R*T)).*eta).*(1 - (F./(2*R*T)).*U);
           
        end
        
        function state = updateSEImassCons(model, state, state0, dt)

            s = model.SEIstochiometricCoeffcient;
            V = model.SEImolarVolume;

            L0 = state0.SEIlength;
            
            L  = state.SEIlength;
            N  = state.SEIflux;

            state.SEImassCons = s/V.*(L - L0)./dt - N;
            
        end
        
        function state = updateSEIvoltageDropEquation(model, state)

            k = model.SEIionicConductivity;

            U = state.SEIvoltageDrop;
            L = state.SEIlength;
            R = state.R;

            state.SEIvoltageDropEquation = U - R.*L.*k;
            
        end

        function state = updateEta(model, state)

            phiElyte = state.phiElectrolyte;
            phiElde  = state.phiElectrode;
            OCP      = state.OCP;
            SEIvoltageDrop     = state.SEIvoltageDrop;

            state.eta = (phiElde - phiElyte - OCP - SEIvoltageDrop);
            
        end

    
    end
end

