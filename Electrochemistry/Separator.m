classdef Separator < PhysicalModel
    
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
            
            model = model@PhysicalModel([]);
            
            fdnames = {'G', ...
                       'thickness', ...
                       'porosity' , ...
                       'rp'       , ...
                       'Gurley'};
            model = dispatchParams(model, paramobj, fdnames);
            
            model.volumeFraction = 1 - model.porosity;

        end
    end
    
end

