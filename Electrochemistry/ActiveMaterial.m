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

        thermalConductivity % Intrinsic Thermal conductivity of the active component
        heatCapacity        % Intrinsic Heat capacity of the active component

        EffectiveThermalConductivity % Effective Thermal Conductivity of the active component
        EffectiveHeatCapacity % Effective Heat Capacity of the active component
        
        electricalConductivity % Electrical conductivity
        
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
            
            paramobj.SolidDiffusion.volumetricSurfaceArea = paramobj.Interface.volumetricSurfaceArea;
            model.SolidDiffusion = SolidDiffusionModel(paramobj.SolidDiffusion);
            
            % setup volumeFraction, porosity, thickness
            nc = model.G.cells.num;
            volumeFraction = model.Interface.volumeFraction*ones(nc, 1);
            model.volumeFraction = volumeFraction;
            model.porosity = 1 - model.volumeFraction;
            model.thickness = 10e-6;
            
            % setup effective electrical conductivity using Bruggeman approximation 
            model.EffectiveElectricalConductivity = model.electricalConductivity.*volumeFraction.^1.5;
            
            % setup effective thermal conductivity            
            model.EffectiveThermalConductivity = model.thermalConductivity.*volumeFraction.^1.5;
            model.EffectiveHeatCapacity = model.heatCapacity.*volumeFraction;
        
            model = model.setupVarPropNames();
                        
        end

        function model = setupVarPropNames(model)

            itf = 'Interface';
            sd = 'SolidDiffusion';
            
            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)
            
            model = model.registerSubModels({'Interface'});
            model = model.registerSubModels({'SolidDiffusion'});
            
            varnames =  {'jCoupling'};
            model = model.registerVarNames(varnames);

            fn = @ActiveMaterial.updatejBcSource;
            model = model.registerPropFunction({'jBcSource', fn, {'jCoupling'}});            

            fn = @ActiveMaterial.updateCurrentSource;
            model = model.registerPropFunction({'eSource', fn, {{itf, 'R'}}});

            fn = @ActiveMaterial.updateChargeCarrier;
            model = model.registerPropFunction({'c', fn, {sd, 'c'}});

            fn = @ActiveMaterial.updatePhi;
            model = model.registerPropFunction({{itf, 'phiElectrode'}, fn, {'phi'}});
            
            fn = @ActiveMaterial.dispatchTemperature;
            model = model.registerPropFunction({{itf, 'T'}, fn, {'T'}});
            model = model.registerPropFunction({{sd, 'T'}, fn, {'T'}});

            fn = @ActiveMaterial.updateSurfaceConcentration;
            model = model.registerPropFunction({{sd, 'cSurface'}, fn, {{sd, 'c'}}});
            model = model.registerPropFunction({{itf, 'cElectrodeSurface'}, fn, {{sd, 'cSurface'}}});
            
            fn = @ActiveMaterial.dispatchSolidRate;
            model = model.registerPropFunction({{sd, 'Rsolid'}, fn, {'R'}});

            
        end
        
        function state = dispatchSolidRate(model, state)
            
            state.SolidDiffusion.Rsolid = state.Interface.R;
            
        end
        
        function state = updateSurfaceConcentration(model, state)

            sd = 'SolidDiffusion';
            itf = 'Interface';
            
            state.(sd) = model.(sd).updateSurfaceConcentration(state.(sd));
            state.(itf).cElectrodeSurface = state.(sd).cSurface;
            
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
