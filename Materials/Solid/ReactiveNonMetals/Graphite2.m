classdef Graphite2 < ActiveElectroChemicalComponent
    
    properties
    end
    
    methods
        
        % Calculate solid diffusion coefficient, [m^2 s^-1]
        D = model.Li.D0 .* exp(-model.Li.EaD./model.con.R*(1./T - 1/refT));
        
    end
    
end
