classdef celgard2500 < FvModel
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
        G       % Gurley number,    [s]
        
        % Mesh properties
        X
        Xb
        N
        dombin
     
        
    end
    
    methods
        function obj = celgard2500(compname, dims, sizes)
            
            obj = obj@FvModel(compname);
            
            obj.t       = 10e-6;
            obj.void    = 0.55;
            obj.eps     = 1 - obj.void;
            obj.rp      = 0.064e-6 ./ 2;
            obj.G       = 200;
            
            obj.Grid = cartGrid(dims, sizes);
            
            obj.varnames = {'phi', 'Li'};
            nc = obj.Grid.cells.num;
            obj.varsizes = [nc, nc];
        end

    end
end

