classdef celgard2500 < ComponentModel
    %UNTITLED10 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        % Physical constants
        con = physicalConstants()
        
        % Physicochemical properties
        t       % Thickness,        [m]
        eps     % Volume fraction,  [-]
        void    % Porosity,         [-]
        rp      % Pore radius,      [m]
        Gurley  % Gurley number,    [s]
        
    end
    
    methods
        function model = celgard2500(name, G, cells)
            
            model = model@ComponentModel(name);
            
            model.t      = 10e-6;
            model.void   = 0.55;
            model.eps    = 1 - model.void;
            model.rp     = 0.064e-6 ./ 2;
            model.Gurley = 200;
            
            model.G = genSubGrid(G, cells);

        end
        
        function state = initializeState(model, state)
            warning('to be implemented');
        end

    end
end

