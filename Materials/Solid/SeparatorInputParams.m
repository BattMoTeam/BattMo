classdef SeparatorInputParams < ComponentInputParams
%
% Input class for :class:`Separator <Electrochemistry.Separator>`
%    
    properties
        
        thickness
        porosity
        rp
        Gurley
        
        thermalConductivity
        heatCapacity

        density % [kg m^-3]
    end
    
    methods

        function paramobj = SeparatorInputParams(jsonstruct)
            paramobj = paramobj@ComponentInputParams(jsonstruct);
        end
        
    end
    
    
end
