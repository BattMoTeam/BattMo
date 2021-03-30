classdef celgard2500 < PhysicalModel
    %UNTITLED10 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        % Physical constants
        con = PhysicalConstants()
        
        % Physicochemical properties
        thickness      % Thickness,        [m]
        volumeFraction % Volume fraction,  [-]
        porosity       % Porosity,         [-]
        rp             % Pore radius,      [m]
        Gurley         % Gurley number,    [s]
        
    end
    
    methods
        function model = celgard2500(G, cells)
            
            model = model@PhysicalModel(G);
            
            model.thickness      = 10e-6;
            model.porosity       = 0.55;
            model.volumeFraction = 1 - model.porosity;
            model.rp             = 0.064e-6 ./ 2;
            model.Gurley         = 200;
            
            model.G = genSubGrid(G, cells);

        end
        
        function state = initializeState(model, state)
            warning('to be implemented');
        end

    end
end

