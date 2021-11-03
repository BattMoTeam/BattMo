classdef ElectrodeActiveComponent < ElectroChemicalComponent
    
    properties
        
        ActiveMaterial % instance of :class:`hello <Electrochemistry.Electrodes.ActiveMaterial>`
        
        volumeFraction % Volume fraction
        porosity       % Porosity
        thickness      % Thickness

        InterDiffusionCoefficient % Inter particle diffusion coefficient parameter (diffusion between the particles)
        
        thermalConductivity % Intrinsic Thermal conductivity
        heatCapacity        % Intrinsic Heat capacity

        EffectiveThermalConductivity % Effective Thermal Conductivity
        EffectiveHeatCapacity % Effective Heat Capacity
        
        electricalConductivity
    end
    
    methods
        
        function model = ElectrodeActiveComponent(paramobj)
        %
        % ``paramobj`` is instance of :class:`ElectrodeActiveComponentInputParams <Electrochemistry.ElectrodeActiveComponentInputParams>`
        %    
            model = model@ElectroChemicalComponent(paramobj);
            
            fdnames = {'thermalConductivity', ...
                       'electricalConductivity', ...
                       'heatCapacity', ...
                       'InterDiffusionCoefficient'};
            model = dispatchParams(model, paramobj, fdnames);
            
            % Setup ActiveMaterial component
            paramobj.ActiveMaterial.G = model.G;
            
            model.ActiveMaterial = ActiveMaterial(paramobj.ActiveMaterial);

            % setup volumeFraction, porosity, thickness
            nc = model.G.cells.num;
            volumeFraction = model.ActiveMaterial.volumeFraction*ones(nc, 1);
            model.volumeFraction = volumeFraction;
            model.porosity = 1 - model.volumeFraction;
            model.thickness = 10e-6;
            
            % setup effective electrical conductivity using Bruggeman approximation 
            model.EffectiveElectricalConductivity = model.electricalConductivity.*volumeFraction.^1.5;
            
            % setup effective diffusion coefficient (inter-particle diffusion)
            model.EffectiveDiffusionCoefficient = model.InterDiffusionCoefficient.*volumeFraction.^1.5;
            
            % setup effective thermal conductivity            
            model.EffectiveThermalConductivity = model.thermalConductivity.*volumeFraction.^1.5;
            model.EffectiveHeatCapacity = model.heatCapacity.*volumeFraction;
            
        end

        function state = updateIonAndCurrentSource(model, state)
            
            F = model.ActiveMaterial.constants.F;
            vols = model.G.cells.volumes;
            n = model.ActiveMaterial.n;
            R = state.ActiveMaterial.R;
            
            state.eSource = - vols.*R*n*F; % C/second
            state.massSource = - vols.*R; % mol/second
            
        end
        
        function state = updateChargeCarrier(model, state)
            state.ActiveMaterial.cElectrodeAveraged = state.c; 
        end
        
        function state = updatePhi(model, state)
            state.ActiveMaterial.phiElectrode = state.phi;
        end         
        
        function state = updateTemperature(model, state)
            state.ActiveMaterial.T = state.T;
        end

        function state = updatejBcSource(model, state)
            state.jBcSource = state.jCoupling;
        end
        
    end
    
end

