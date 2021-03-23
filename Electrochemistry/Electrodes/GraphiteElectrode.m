classdef GraphiteElectrode < Electrode
    
    % GRAPHITE ELECTRODE Electrode class with graphite active material.
    
    properties

        SEI % Solid-electrolyte interphase (SEI) object
        
    end
    
    methods
        
        function model = GraphiteElectrode(G, cells)
            
            model = model@Electrode();
            
            % setup grid
            G = genSubGrid(G, cells);
            model.G = G;
            
            % setup discrete differential operators
            model.operators = localSetupOperators(G);
            
            % setup graphite active material
            ActiveMaterial = Graphite();
            model.ActiveMaterial = ActiveMaterial;
            
            % model.bin  = ptfe();
            % model.sei  = seiAM();
            nc = G.cells.num;
            volumeFraction = ActiveMaterial.volumeFraction*ones(nc, 1);
            model.volumeFraction = volumeFraction;
            model.porosity = 1 - model.volumeFraction;
            model.thickness = 10e-6;
            
            % setup effective electronic conductivity
            econd = ActiveMaterial.electronicConductivity;
            model.EffectiveElectronicConductivity = econd .* volumeFraction.^1.5;
            
        end
        
        function state = updateReactionRate(model, state)
        % Abbreviation used in this function
        % am : ActiveMaterial
            
            am = model.ActiveMaterial;
            
            T = state.T;
            phiElyte = state.phiElectrolyte;
            
            phi = state.ActiveMaterial.phi;
            OCP = state.ActiveMaterial.OCP;
            k = state.ActiveMaterial.k;
            
            eta = (phi - phiElyte - OCP);
            F = model.constants.F;
            R = am.volumetricSurfaceArea.*ButlerVolmerEquation(k.*model.constants.F, 0.5, 1, eta, T)/F;
             
            state.eta = eta;
            state.R = R;
            
        end

        function state = updateIonFlux(model, state)
            state = assembleLithiumFlux(model, state);
        end 
        
    end
end

                        
