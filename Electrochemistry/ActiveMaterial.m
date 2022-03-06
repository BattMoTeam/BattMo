classdef ActiveMaterial < ElectronicComponent
    
    properties
        
        %
        % instance of :class:`Interface <Electrochemistry.Electrodes.Interface>`
        %
        Interface 
        SolidDiffusion
        
        volumeFraction % Volume fraction
        porosity       % Porosity
        thickness      % Thickness / [m]
        
        InterDiffusionCoefficient % Inter particle diffusion coefficient parameter (diffusion between the particles)
        EffectiveDiffusionCoefficient % 
        
        thermalConductivity % Intrinsic Thermal conductivity of the active component
        heatCapacity        % Intrinsic Heat capacity of the active component

        EffectiveThermalConductivity % Effective Thermal Conductivity of the active component
        EffectiveHeatCapacity % Effective Heat Capacity of the active component
        
        electricalConductivity % Electrical conductivity

        useSimplifiedDiffusionModel
    end
    
    methods
        
        function model = ActiveMaterial(paramobj)
        %
        % ``paramobj`` is instance of :class:`ActiveMaterialInputParams <Electrochemistry.ActiveMaterialInputParams>`
        %    
            model = model@ElectronicComponent(paramobj);

            fdnames = {'thermalConductivity'   , ...
                       'electricalConductivity', ...
                       'heatCapacity'};
            
            model = dispatchParams(model, paramobj, fdnames);
            
            % Setup Interface component
            paramobj.Interface.G = model.G;
            
            model.Interface = Interface(paramobj.Interface);
            
            paramobj.SolidDiffusion.volumetricSurfaceArea = model.Interface.volumetricSurfaceArea;
            
            if paramobj.SolidDiffusion.useSimplifiedDiffusionModel
                model.SolidDiffusion = SimplifiedSolidDiffusionModel(paramobj.SolidDiffusion);
                model.useSimplifiedDiffusionModel = true;
            else
                model.SolidDiffusion = SolidDiffusionModel(paramobj.SolidDiffusion);
                model.useSimplifiedDiffusionModel = false;
            end
            
            if model.useSimplifiedDiffusionModel
                % only used in SimplifiedSolidDiffusionModel (for now)
                model.InterDiffusionCoefficient = paramobj.InterDiffusionCoefficient;
            end
            
            % setup volumeFraction, porosity, thickness
            nc = model.G.cells.num;
            volumeFraction = model.Interface.volumeFraction*ones(nc, 1);
            model.volumeFraction = volumeFraction;
            model.porosity = 1 - model.volumeFraction;
            model.thickness = 10e-6;
            
            % setup effective electrical conductivity using Bruggeman approximation 
            model.EffectiveElectricalConductivity = model.electricalConductivity.*volumeFraction.^1.5;
            
            if model.useSimplifiedDiffusionModel            
                model.EffectiveDiffusionCoefficient = model.InterDiffusionCoefficient.*volumeFraction.^1.5;
            end
            
            % setup effective thermal conductivity            
            model.EffectiveThermalConductivity = model.thermalConductivity.*volumeFraction.^1.5;
            model.EffectiveHeatCapacity = model.heatCapacity.*volumeFraction;
        
            model = model.setupVarPropNames();
                        
        end

        function model = setupVarPropNames(model)

            itf = 'Interface';
            sd  = 'SolidDiffusion';
            
            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)
            
            model = model.registerSubModels({'Interface'});
            model = model.registerSubModels({'SolidDiffusion'});
            
            varnames =  {'jCoupling'};
            if  model.useSimplifiedDiffusionModel
                varnames{end + 1} = {'c'};
                varnames{end + 1} = {'massCons'};
                varnames{end + 1} = {'massAccum'};
                varnames{end + 1} = {'massSource'};
            end
            model = model.registerVarNames(varnames);

            fn = @ActiveMaterial.updatejBcSource;
            model = model.registerPropFunction({'jBcSource', fn, {'jCoupling'}});            

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
        
        function state = dispatchRate(model, state)
            
            state.SolidDiffusion.R = state.Interface.R;
            
        end
        
        function state = updateConcentrations(model, state)

            sd  = 'SolidDiffusion';
            itf = 'Interface';
            
            if model.useSimplifiedDiffusionModel
                state.(sd).cAverage = state.c;
                state.(itf).cElectrodeSurface = state.(sd).cSurface;
            else
                state.(sd) = model.(sd).updateSurfaceConcentration(state.(sd));
                state.(itf).cElectrodeSurface = state.(sd).cSurface;
            end
            
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

        function state = addSOC(model, state)

            itf = model.Interface; 
            c = state.c; 
            
            theta = c/itf.Li.cmax; 
            m = (1 ./ (itf.theta100 - itf.theta0)); 
            b = - m.*itf.theta0; 
            
            state.SOC = theta*m + b;
            state.theta = theta;
            
        end
        
    end
    
end





%{
Copyright 2009-2021 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BatMo

BatMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BatMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BatMo.  If not, see <http://www.gnu.org/licenses/>.
%}
