classdef ElectrodeActiveComponent < ElectroChemicalComponent
    
    properties
        
        ActiveMaterial
        
        volumeFraction
        porosity
        thickness
                    
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
                       'heatCapacity'};
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
            % setup effective thermal conductivity            
            model.EffectiveThermalConductivity = model.thermalConductivity.*volumeFraction.^1.5;
            model.EffectiveHeatCapacity = model.heatCapacity.*volumeFraction;
            
        end

        function state = updateIonAndCurrentSource(model, state)
            
            ccSourceName = model.chargeCarrierSourceName;
            F = model.ActiveMaterial.constants.F;
            vols = model.G.cells.volumes;
            
            R = state.ActiveMaterial.R;
            
            state.eSource = - vols.*R;
            state.(ccSourceName) = - vols.*R/F;
            
        end
        
        function state = updateDiffusionCoefficient(model, state)
            
            D = state.ActiveMaterial.D;
            state.D = D .* model.volumeFraction .^1.5;
            
        end
        
        function state = updateChargeCarrier(model, state)
            state.ActiveMaterial.cElectrode = state.c; 
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

