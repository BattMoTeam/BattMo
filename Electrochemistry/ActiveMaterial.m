classdef ActiveMaterial < ElectronicComponent
    
    properties
        
        %
        % instance of :class:`Interface <Electrochemistry.Electrodes.Interface>`
        %
        Interface 
        
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
            
            % defines shortcuts
            am = 'ActiveMaterial';
            cc  = 'CurrentCollector';
            am  = 'Interface';
            sd  = 'SolidDiffusion';
            
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
        
            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)
            
            model = model.registerSubModels({'Interface'});

            varnames =  {'jCoupling'};
            model = model.registerVarNames(varnames);

            fn = @ActiveMaterial.updatejBcSource;
            model = model.registerPropFunction({'jBcSource', fn, {'jCoupling'}});            

            fn = @ActiveMaterial.updateCurrentSource;
            model = model.registerPropFunction({'eSource', fn, {{am, 'R'}}});

            fn = @ActiveMaterial.updateChargeCarrier;
            model = model.registerPropFunction({'c', fn, {am, sd, 'c'}});

            fn = @ActiveMaterial.updatePhi;
            model = model.registerPropFunction({{am, 'phiElectrode'}, fn, {'phi'}});
            
            fn = @ActiveMaterial.updateTemperature;
            model = model.registerPropFunction({{am, 'T'}, fn, {'T'}});
            model = model.registerPropFunction({{am, sd, 'T'}, fn, {'T'}});
            
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
        
        function state = updateTemperature(model, state)
            state.Interface.T = state.T;
            state.Interface.SolidDiffusion.T = state.T;
        end

        function state = updatejBcSource(model, state)
            state.jBcSource = state.jCoupling;
            state.jFaceBc = state.jFaceCoupling;
        end

        function state = addSOC(model, state)

            am = model.Interface; 
            c = state.c; 
            
            theta = c/am.Li.cmax; 
            m = (1 ./ (am.theta100 - am.theta0)); 
            b = - m.*am.theta0; 
            
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
