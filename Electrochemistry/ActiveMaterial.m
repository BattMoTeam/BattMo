classdef ActiveMaterial < ElectronicComponent
    
    properties
        
        %
        % instance of :class:`Interface <Electrochemistry.Electrodes.Interface>`
        %


        Interface

        SolidDiffusion        

        porosity                      % porosity
        volumeFraction                % Volume fraction
        electricalConductivity        % Electrical conductivite
        InterDiffusionCoefficient     % Inter particle diffusion coefficient parameter (diffusion between the particles)
        thermalConductivity           % Intrinsic Thermal conductivity of the active component
        heatCapacity                  % Intrinsic Heat capacity of the active component

        EffectiveDiffusionCoefficient % 

        EffectiveThermalConductivity  % Effective Thermal Conductivity of the active component
        EffectiveHeatCapacity         % Effective Heat Capacity of the active component

        useSimplifiedDiffusionModel
        
        externalCouplingTerm          % only used in no current collector

    end
    
    methods
        
        function model = ActiveMaterial(paramobj)
        %
        % ``paramobj`` is instance of :class:`ActiveMaterialInputParams <Electrochemistry.ActiveMaterialInputParams>`
        %    
            model = model@ElectronicComponent(paramobj);

            fdnames = {'thermalConductivity'   , ...
                       'electricalConductivity', ...
                       'heatCapacity', ...
                       'externalCouplingTerm', ...
                       'useSimplifiedDiffusionModel'};
            
            model = dispatchParams(model, paramobj, fdnames);
            
            % Setup Interface component
            paramobj.Interface.G = model.G;
            
            model.Interface = Interface(paramobj.Interface);
            
            paramobj.SolidDiffusion.volumetricSurfaceArea = model.Interface.volumetricSurfaceArea;
            
            if paramobj.useSimplifiedDiffusionModel
                model.SolidDiffusion = SimplifiedSolidDiffusionModel(paramobj.SolidDiffusion);
            else
                paramobj.SolidDiffusion.np = model.G.cells.num;
                model.SolidDiffusion = FullSolidDiffusionModel(paramobj.SolidDiffusion);
            end
            
            if model.useSimplifiedDiffusionModel
                % only used in SimplifiedSolidDiffusionModel (for now)
                model.InterDiffusionCoefficient = paramobj.InterDiffusionCoefficient;
            end
            
            nc = model.G.cells.num;

            volumeFraction = model.Interface.volumeFraction*ones(nc, 1);
            model.porosity = 1 - volumeFraction;

            model = model.setupDependentProperties();
            
        end
        
        function model = setupDependentProperties(model)           

            model.volumeFraction = 1 - model.porosity;
            volumeFraction = model.volumeFraction;
            % setup effective electrical conductivity using Bruggeman approximation 
            model.EffectiveElectricalConductivity = model.electricalConductivity.*volumeFraction.^1.5;
            
            if model.useSimplifiedDiffusionModel            
                model.EffectiveDiffusionCoefficient = model.InterDiffusionCoefficient.*volumeFraction.^1.5;
            end
            
            % setup effective thermal conductivity            
            model.EffectiveThermalConductivity = model.thermalConductivity.*volumeFraction.^1.5;
            model.EffectiveHeatCapacity = model.heatCapacity.*volumeFraction;
        end

        function model = registerVarAndPropfuncNames(model)

            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@ElectronicComponent(model);
            
            itf = 'Interface';
            sd  = 'SolidDiffusion';
            
            varnames =  {'jCoupling'};
            if  model.useSimplifiedDiffusionModel
                varnames{end + 1} = {'c'};
                varnames{end + 1} = {'massCons'};
                varnames{end + 1} = {'massAccum'};
                varnames{end + 1} = {'massSource'};
            end
            model = model.registerVarNames(varnames);

            isRoot = true;
            if isRoot
                fn = @ActiveMaterial.updateStandalonejBcSource;
                model = model.registerPropFunction({'jBcSource', fn, {'controlCurrentSource'}});            
            else
                fn = @ActiveMaterial.updatejBcSource;
                model = model.registerPropFunction({'jBcSource', fn, {'jCoupling', 'jFaceCoupling'}});            
            end
            
            fn = @ActiveMaterial.updateCurrentSource;
            model = model.registerPropFunction({'eSource', fn, {{itf, 'R'}}});
            
            fn = @ActiveMaterial.updatePhi;
            model = model.registerPropFunction({{itf, 'phiElectrode'}, fn, {'phi'}});
            
            fn = @ActiveMaterial.dispatchTemperature;
            model = model.registerPropFunction({{itf, 'T'}, fn, {'T'}});
            model = model.registerPropFunction({{sd, 'T'}, fn, {'T'}});

            fn = @ActiveMaterial.updateConcentrations;
            if model.useSimplifiedDiffusionModel
                model = model.registerPropFunction({{sd, 'cAverage'}, fn, {'c'}});
                model = model.registerPropFunction({{itf, 'cElectrodeSurface'}, fn, {{sd, 'cSurface'}}});
            else
                model = model.registerPropFunction({{sd, 'cSurface'}, fn, {{sd, 'c'}}});
                model = model.registerPropFunction({{itf, 'cElectrodeSurface'}, fn, {{sd, 'cSurface'}}});
            end

            if  model.useSimplifiedDiffusionModel
                fn = @ActiveMaterial.updateMassFlux;
                model = model.registerPropFunction({'massFlux', fn, {'c'}});
                fn = @ActiveMaterial.updateMassSource;
                model = model.registerPropFunction({'massSource', fn, {{itf, 'R'}}});
                fn = @ActiveMaterial.updateMassConservation;
                model = model.registerPropFunction({'massCons', fn, {'massAccum', 'massSource'}});
            end
            
            fn = @ActiveMaterial.dispatchRate;
            model = model.registerPropFunction({{sd, 'R'}, fn, {{itf, 'R'}}});
            
            
        end
        
        
        function [problem, state] = getEquations(model, state0, state,dt, drivingForces, varargin)
            
            sd = 'SolidDiffusion';
            itf = 'Interface';          
            
            time = state0.time + dt;
            state = model.initStateAD(state);
            
            state = model.dispatchTemperature(state);
            state = model.updateConcentrations(state);

            %% Assemble reaction rates
            state.(itf) = model.(itf).updateOCP(state.(itf));
            state       = model.updatePhi(state);
            state.(itf) = model.(itf).updateReactionRateCoefficient(state.(itf));
            state.(itf) = model.(itf).updateReactionRate(state.(itf));

            %% Setup charge conservation equation
            
            state = model.updateControl(state, drivingForces, dt);
            state = model.updateCurrent(state);
            state = model.updateCurrentSource(state);
            state = model.updateStandalonejBcSource(state);
            state = model.updateChargeConservation(state);
            
            %% Setup mass conservation equation
            state.(sd) = model.(sd).updateDiffusionCoefficient(state.(sd));
            state.(sd) = model.(sd).updateFlux(state.(sd));
            state      = model.dispatchRate(state);
            state.(sd) = model.(sd).updateMassSource(state.(sd));
            state.(sd) = model.(sd).updateAccumTerm(state.(sd), state0.(sd), dt);
            state.(sd) = model.(sd).updateMassConservation(state.(sd));

            state.(sd) = model.(sd).assembleSolidDiffusionEquation(state.(sd));
            
            %% Setup equations and add some scaling
            %% FIXME: get robust scalings
            massConsScaling = model.(sd).constants.F; % Faraday constant
            eqs = {};
            eqs{end + 1} = state.chargeCons;
            eqs{end + 1} = 1e18*state.(sd).massCons;
            eqs{end + 1} = 1e5*state.(sd).solidDiffusionEq.*massConsScaling.*model.(itf).G.cells.volumes/dt;
            
            names = {'chargeCons', ...
                     'massCons', ...
                     'solidDiffusionEq'};
            
            types = {'cell', 'cell', 'cell'};

            primaryVars = model.getPrimaryVariables();
            
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

        end

        %% Functions for Newton API
        
        function primaryvarnames = getPrimaryVariables(model)
            
            sd = 'SolidDiffusion';
            
            primaryvarnames = {{sd, 'c'}, ...
                               {sd, 'cSurface'}, ...
                               {'phi'}};
            
        end
        
        function forces = getValidDrivingForces(model)

            forces = getValidDrivingForces@PhysicalModel(model);
            forces.src = [];
        end

        function state = updateControl(model, state, drivingForces, dt)
            
            G = model.G;
            coef = G.cells.volumes;
            coef = coef./(sum(coef));
            
            state.controlCurrentSource = drivingForces.src.*coef;
            
        end
        
        
        function cleanState = addStaticVariables(model, cleanState, state, state0)
            
            cleanState = addStaticVariables@BaseModel(model, cleanState, state);
            
            itf = 'Interface';
            
            cleanState.T = state.T;
            cleanState.(itf).cElectrolyte = state.(itf).cElectrolyte;
            cleanState.(itf).phiElectrolyte = state.(itf).phiElectrolyte;
            
        end

        
        function [state, report] = updateState(model, state, problem, dx, drivingForces)

            [state, report] = updateState@BaseModel(model, state, problem, dx, drivingForces);
            
        end
        
        function [model, state] = prepareTimestep(model, state, state0, dt, drivingForces)
            
            [model, state] = prepareTimestep@BaseModel(model, state, state0, dt, drivingForces);
            
        end
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            
            % [state, report] = updateAfterConvergence@BaseModel(model, state0, state, dt, drivingForces);
            report = [];
        end
         
        function model = validateModel(model, varargin)
            
        end

        %% assembly functions use in this model
         
        function state = dispatchRate(model, state)
            
            state.SolidDiffusion.R = state.Interface.R;
            
        end
        
        function state = updateConcentrations(model, state)

            sd  = 'SolidDiffusion';
            itf = 'Interface';
            
            if model.useSimplifiedDiffusionModel
                state.(sd).cAverage = state.c;
            end
            
            state.(itf).cElectrodeSurface = state.(sd).cSurface;
            
        end

        function state = updateMassFlux(model, state)
        % used when useSimplifiedDiffusionModel is true

            D = model.EffectiveDiffusionCoefficient;
            
            c = state.c;

            massflux = assembleFlux(model, c, D);
            
            state.massFlux = massflux;

        end
            
        function state = assembleAccumTerm(model, state, state0, dt)
        % used when useSimplifiedDiffusionModel is true
            
            volumeFraction = model.volumeFraction;
            vols = model.G.cells.volumes;
            
            c  = state.c;
            c0 = state0.c;
            
            state.massAccum = vols.*volumeFraction.*(c - c0)/dt;
            
        end

        function state = updateMassSource(model, state)
        % used when useSimplifiedDiffusionModel is true
            
            vols = model.G.cells.volumes;
            
            R = state.Interface.R;
            
            state.massSource = - R.*vols;
            
        end
        
        function state = updateMassConservation(model, state)
        % used when useSimplifiedDiffusionModel is true
        % Here, we have no flux (for the moment)
            
            state.massCons = state.massAccum - state.massSource;
            
        end
        
        function state = updateStandalonejBcSource(model, state)
            
            state.jBcSource = state.controlCurrentSource;

        end

        function state = updateCurrentSource(model, state)
            
            F = model.Interface.constants.F;
            vols = model.G.cells.volumes;
            n = model.Interface.n;

            R = state.Interface.R;
            
            state.eSource = - vols.*R*n*F; % C/second
            
        end
        
        function state = updatePhi(model, state)
            state.Interface.phiElectrode = state.phi;
        end         
        
        function state = dispatchTemperature(model, state)
            state.Interface.T = state.T;
            state.SolidDiffusion.T = state.T;
        end

        function state = updatejBcSource(model, state)
            state.jBcSource = state.jCoupling;
            state.jFaceBc = state.jFaceCoupling;
        end
        
        function state = updatejBcSourceNoCurrentCollector(model, state)
            state.jBcSource = state.jExternal;
            state.jFaceBc   = state.jFaceExternal;
        end

        function state = addSOC(model, state)

            itf = model.Interface; 
            c = state.c; 
            
            theta = c/itf.cmax; 
            m = (1 ./ (itf.theta100 - itf.theta0)); 
            b = - m.*itf.theta0; 
            
            state.SOC = theta*m + b;
            state.theta = theta;
            
        end
        
        
    end
    
end


%{
Copyright 2021-2022 SINTEF Industry, Sustainable Energy Technology
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
