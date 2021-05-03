classdef Separator < ThermalComponent.m
    
    properties
        
        % Physicochemical properties
        thickness      % Thickness,        [m]
        porosity       % Porosity,         [-]
        rp             % Pore radius,      [m]
        Gurley         % Gurley number,    [s]
        volumeFraction % Volume fraction,  [-]
        
    end
    
    methods
        function model = Separator(paramobj)
            
            model = model@ThermalComponent(paramobj);
            
            fdnames = {'thermalConductivity', ...
                       'thickness', ...
                       'porosity' , ...
                       'rp'       , ...
                       'Gurley'};
            model = dispatchParams(model, paramobj, fdnames);
            
            model.volumeFraction = 1 - model.porosity;
            model.EffectiveThermalConductivity = model.volumeFraction.*model.thermalConductivity;

        end
    end
    
end

