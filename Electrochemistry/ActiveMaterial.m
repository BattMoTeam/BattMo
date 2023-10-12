classdef ActiveMaterial < BaseModel
    
    properties
        
        %
        % instance of :class:`Interface <Electrochemistry.Electrodes.Interface>`
        %

        %% Sub-Models
        
        Interface
        SolidDiffusion        

        %% Input parameters

        % Standard parameters

        electronicConductivity % the electronic conductivity of the material (symbol: sigma)
        density                % the mass density of the material (symbol: rho)
        massFraction           % the ratio of the mass of the material to the total mass of the phase or mixture (symbol: gamma)
        
        thermalConductivity    % the intrinsic Thermal conductivity of the active component
        specificHeatCapacity   % Specific Heat capacity of the active component

        diffusionModelType     % either 'full' or 'simple'
        
        % Advanded parameters
        
        standAlone % Set to true if model is used as main model for development purpose (standard use is as a sub-model)

        % Coupling parameters
        
        externalCouplingTerm % structure to describe external coupling (used in absence of current collector)


        % Others
        
        % properties needed when model is run in stand-alone mode
        
        funcCallList
        primaryVarNames
        equationVarNames
        equationNames
        equationTypes
        

        
    end
    
    methods
        
        function model = ActiveMaterial(paramobj)
        %
        % ``paramobj`` is instance of :class:`ActiveMaterialInputParams <Electrochemistry.ActiveMaterialInputParams>`
        %
            model = model@BaseModel();
            
            fdnames = {'electronicConductivity', ... 
                       'density'               , ...                
                       'massFraction'          , ...
                       'thermalConductivity'   , ...    
                       'specificHeatCapacity'  , ...   
                       'volumeFraction'        , ... 
                       'externalCouplingTerm'  , ...
                       'diffusionModelType'    , ...
                       'standAlone'};

            model = dispatchParams(model, paramobj, fdnames);

            model.Interface = Interface(paramobj.Interface);

            diffusionModelType = model.diffusionModelType;

            switch model.diffusionModelType
              case 'simple'
                model.SolidDiffusion = SimplifiedSolidDiffusionModel(paramobj.SolidDiffusion);
              case 'full'
                model.SolidDiffusion = FullSolidDiffusionModel(paramobj.SolidDiffusion);
              otherwise
                error('Unknown diffusionModelType %s', diffusionModelType);
            end

            if model.standAlone

                model = model.setupStandAloneModel();
                
            end
            
            
        end
        
        function model = registerVarAndPropfuncNames(model)

            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@BaseModel(model);
            
            itf = 'Interface';
            sd  = 'SolidDiffusion';

            varnames = {'T'};
            model = model.registerVarNames(varnames);

            fn = @ActiveMaterial.dispatchTemperature;
            model = model.registerPropFunction({{sd, 'T'}, fn, {'T'}});
            model = model.registerPropFunction({{itf, 'T'}, fn, {'T'}});
            
            if model.standAlone

                varnames = {};

                % Current density
                varnames{end + 1} = 'j';
                % Ootential in Electrode
                varnames{end + 1} = 'phi';
                % Charge Conservation equation
                varnames{end + 1} = 'chargeCons';
                
                model = model.registerVarNames(varnames);

                varnames = {{itf, 'dUdT'}, ...
                            'jCoupling', ...
                            'jExternal'};
                model = model.removeVarNames(varnames);

                varnames = {'T', ...
                            {itf, 'cElectrolyte'},... 
                            {itf, 'phiElectrolyte'}};
                model = model.registerStaticVarNames(varnames);
                
            end

            fn = @ActiveMaterial.updateRvol;
            model = model.registerPropFunction({{sd, 'Rvol'}, fn, {{itf, 'R'}}});
            
            fn = @ActiveMaterial.updateConcentrations;
            model = model.registerPropFunction({{itf, 'cElectrodeSurface'}, fn, {{sd, 'cSurface'}}});
            
            if model.standAlone
                
                fn = @ActiveMaterial.updateControl;
                fn = {fn, @(propfunction) PropFunction.drivingForceFuncCallSetupFn(propfunction)};
                model = model.registerPropFunction({'j', fn, {}});
                
                fn = @ActiveMaterial.updateChargeCons;
                inputnames = {'j', ...
                              {sd, 'Rvol'}}
                model = model.registerPropFunction({'chargeCons', fn, inputnamse});

            end
            
        end
        

        function model = setupStandAloneModel(model)

            model = model.setupComputationalGraph();

            % cgt = model.computationalGraph();
            
            % model.funcCallList     = cgt.getOrderedFunctionCallList();
            % model.primaryVarNames  = cgt.getPrimaryVariableNames();
            % model.equationVarNames = cgt.getEquationVariableNames();
            
        end

        
        function [problem, state] = getEquations(model, state0, state,dt, drivingForces, varargin)
            
            sd  = 'SolidDiffusion';
            itf = 'Interface';
            
            time = state0.time + dt;
            state = model.initStateAD(state);

            state = updateControl(model, state, drivingForces);
            
            state                = model.updateStandalonejBcSource(state);
            state                = model.updateCurrent(state);
            state.SolidDiffusion = model.SolidDiffusion.updateMassAccum(state.SolidDiffusion, state0.SolidDiffusion, dt);
            state                = model.dispatchTemperature(state);
            state.SolidDiffusion = model.SolidDiffusion.updateDiffusionCoefficient(state.SolidDiffusion);
            state.SolidDiffusion = model.SolidDiffusion.updateFlux(state.SolidDiffusion);
            state                = model.updateConcentrations(state);
            state                = model.updatePhi(state);
            state.Interface      = model.Interface.updateReactionRateCoefficient(state.Interface);
            state.Interface      = model.Interface.updateOCP(state.Interface);
            state.Interface      = model.Interface.updateEta(state.Interface);
            state.Interface      = model.Interface.updateReactionRate(state.Interface);
            state                = model.updateRvol(state);
            state                = model.updateCurrentSource(state);
            state                = model.updateChargeConservation(state);
            state.SolidDiffusion = model.SolidDiffusion.updateMassSource(state.SolidDiffusion);
            state.SolidDiffusion = model.SolidDiffusion.assembleSolidDiffusionEquation(state.SolidDiffusion);
            state.SolidDiffusion = model.SolidDiffusion.updateMassConservation(state.SolidDiffusion);
            
            %% Setup equations and add some scaling
            n     = model.(itf).n; % number of electron transfer (equal to 1 for Lithium)
            F     = model.(sd).constants.F;
            vol   = model.operators.pv;
            rp    = model.(sd).rp;
            vsf   = model.(sd).volumetricSurfaceArea;
            surfp = 4*pi*rp^2;
            
            scalingcoef = (vsf*vol(1)*n*F)/surfp;
            
            eqs = {};
            eqs{end + 1} = state.chargeCons;
            eqs{end + 1} = scalingcoef*state.(sd).massCons;
            eqs{end + 1} = scalingcoef*state.(sd).solidDiffusionEq;
            
            names = {'chargeCons', ...
                     'massCons', ...
                     'solidDiffusionEq'};
            
            types = {'cell', 'cell', 'cell'};

            primaryVars = model.getPrimaryVariableNames();
            
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

        end

        function primaryvarnames = getPrimaryVariableNames(model)
            
            sd = 'SolidDiffusion';
            
            primaryvarnames = {{sd, 'c'}, ...
                               {sd, 'cSurface'}, ...
                               {'phi'}};
            
        end
        
        function forces = getValidDrivingForces(model)

            forces = getValidDrivingForces@PhysicalModel(model);
            forces.src = [];
            
        end

        function state = updateControl(model, state, drivingForces)
            
            G = model.G;
            coef = G.cells.volumes;
            coef = coef./(sum(coef));
            
            state.controlCurrentSource = drivingForces.src.*coef;
            
        end
        
        
        function cleanState = addStaticVariables(model, cleanState, state, state0)
            
            cleanState = addStaticVariables@BaseModel(model, cleanState, state);
            
            itf = 'Interface';
            
            cleanState.T = state.T;
            cleanState.(itf).cElectrolyte   = state.(itf).cElectrolyte;
            cleanState.(itf).phiElectrolyte = state.(itf).phiElectrolyte;
            
            sigma = model.electronicConductivity;
            vf    = model.volumeFraction;
            brugg = model.bruggemanCoefficient;
            
            cleanState.conductivity = sigma*vf.^brugg;
            
        end

        
        function [state, report] = updateState(model, state, problem, dx, drivingForces)

            [state, report] = updateState@BaseModel(model, state, problem, dx, drivingForces);
            
        end
        
        function [model, state] = prepareTimestep(model, state, state0, dt, drivingForces)
            
            [model, state] = prepareTimestep@BaseModel(model, state, state0, dt, drivingForces);
            
        end
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
        % [state, report] = updateAfterConvergence@ElectronicComponent(model, state0, state, dt, drivingForces);

        % by not calling the parent method, we do not clean the state s
            report = [];
            
        end
         
        function model = validateModel(model, varargin)
            
        end

        %% assembly functions use in this model

        function state = updateRvol(model, state)
            
            itf = 'Interface';
            sd  = 'SolidDiffusion';
            
            vsa = model.(itf).volumetricSurfaceArea;
            
            Rvol = vsa.*state.(itf).R;
            
            state.(sd).Rvol = Rvol;
            
        end        
        
        function state = updateConcentrations(model, state)

            sd  = 'SolidDiffusion';
            itf = 'Interface';
            
            state.(itf).cElectrodeSurface = state.(sd).cSurface;
            
        end

        
        function state = updateStandalonejBcSource(model, state)
            
            state.jBcSource = state.controlCurrentSource;

        end

        
        function state = dispatchTemperature(model, state)
            state.Interface.T = state.T;
            state.SolidDiffusion.T = state.T;
        end



        function state = updateAverageConcentration(model, state)

            % shortcut
            sd  = 'SolidDiffusion';

            vf       = model.volumeFraction;
            am_frac  = model.activeMaterialFraction;
            vols     = model.G.cells.volumes;
            
            c = state.(sd).cAverage;

            vols = am_frac*vf.*vols;

            cAverage = sum(c.*vols)/sum(vols);

            state.cAverage = cAverage;
            
        end
        
        
        
    end
    
end


%{
Copyright 2021-2023 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
