classdef ElectrodeActiveComponent < ElectroChemicalComponent
% reference used here
% @article{ref1,
% 	year = 2007,
% 	author = {Qi Zhang and Ralph E. White},
% 	title = {Comparison of approximate solution methods for the solid phase diffusion equation in a porous electrode model},
% 	journal = {Journal of Power Sources}
% }
    
    properties
        
        ActiveMaterial
        
        volumeFraction
        porosity
        thickness
                
        % Inter particle diffusion coefficient parameter (diffusion between the particles)
        InterDiffusionCoefficient
        
        thermalConductivity % intrinsic thermal conductivity value
        heatCapacity % intrinsic heat capacity value

        EffectiveThermalConductivity
        EffectiveHeatCapacity
        
    end
    
    methods
        
        function model = ElectrodeActiveComponent(paramobj)
        % shortcut used here:
        % am = ActiveMaterial
            
            model = model@ElectroChemicalComponent(paramobj);
            
            fdnames = {'thermalConductivity', ...
                       'heatCapacity', ...
                       'InterDiffusionCoefficient'};
            model = dispatchParams(model, paramobj, fdnames);
            
            % Setup ActiveMaterial component
            paramobj.am.G = model.G;
            
            model.ActiveMaterial = setupActiveMaterial(paramobj.am);

            % setup volumeFraction, porosity, thickness
            nc = model.G.cells.num;
            volumeFraction = model.ActiveMaterial.volumeFraction*ones(nc, 1);
            model.volumeFraction = volumeFraction;
            model.porosity = 1 - model.volumeFraction;
            model.thickness = 10e-6;
            
            % setup effective electrical conductivity
            econd = model.ActiveMaterial.electricalConductivity;
            % Bruggeman approximation 
            model.EffectiveElectricalConductivity = econd .* volumeFraction.^1.5;
            
            % setup effective diffusion coefficient (inter-particle diffusion)
            model.EffectiveDiffusionCoefficient = model.InterDiffusionCoefficient.*volumeFraction.^1.5;
            
            % setup effective thermal conductivity            
            model.EffectiveThermalConductivity = model.thermalConductivity.*volumeFraction.^1.5;
            model.EffectiveHeatCapacity = model.heatCapacity.*volumeFraction;
            
        end

        function state = updateIonAndCurrentSource(model, state)
            
            ccSourceName = model.chargeCarrierSourceName;
            F = model.ActiveMaterial.constants.F;
            vols = model.G.cells.volumes;
            n = model.ActiveMaterial.n;
            R = state.ActiveMaterial.R;
            
            state.eSource = - vols.*R*n*F; % C/second
            state.(ccSourceName) = - vols.*R; % mol/second
            
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

