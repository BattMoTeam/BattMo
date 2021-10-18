classdef BinaryElectrolyte < Electrolyte
%
% Binary Electrolyte model
    
    properties

        diffusionConcentrationFittingCoefficients
        diffusionTemperatureFittingCoefficients
        
    end
    
    methods
        
        function model = BinaryElectrolyte(paramobj)
            
            model = model@Electrolyte(paramobj);
            fdnames = {'ionicConductivityFittingCoefficients', 
                       'diffusionConcentrationFittingCoefficients',
                       'diffusionTemperatureFittingCoefficients'};
            model = dispatchParams(model, paramobj, fdnames);
            
        end
        
        function state = updateConcentrations(model, state)
            state.cs{2} = state.cs{1};
        end


        function state = updateDiffusionCoefficient(model, state)
            
            c = state.cs{1};
            T = state.T;
            
            % Calculate diffusion coefficients constant for the diffusion coefficient calcuation
            cnst = model.diffusionConcentrationFittingCoefficients;
            Tgi = model.diffusionTemperatureFittingCoefficients;
            
            % Diffusion coefficient, [m^2 s^-1]
            D = 1e-4 .* 10 .^ ( ( cnst(1,1) + cnst(1,2) ./ ( T - Tgi(1) - Tgi(2) .* c .* 1e-3) + cnst(2, 1) .* ...
                                  c .* 1e-3) );
            
            % calculate the effective diffusion coefficients in porous media
            state.D = D .* model.volumeFraction .^1.5;
        
        end
        
    end
end

