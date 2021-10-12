classdef BinaryElectrolyteInputParams < ElectrolyteInputParams
%
% Input class for :class:`BinaryElectrolyte <Electrochemistry.Electrolyte>`
%    
    properties
        
        ionicConductivityFittingCoefficients
        diffusionConcentrationFittingCoefficients
        diffusionTemperatureFittingCoefficients
        
    end
    
    methods

        function paramobj = BinaryElectrolyteInputParams(jsonstruct);
            paramobj = paramobj@ElectrolyteInputParams(jsonstruct);
        end

    end
    
end
