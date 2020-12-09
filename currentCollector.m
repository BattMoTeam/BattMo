classdef currentCollector < FvModel
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        % Physical constants
        con = physicalConstants();
        
        % Design properties
        t % Thickness,        [m]
        
        % Material properties
        am % Active material object
        
        % State properties
        j % Current density,      [A/m2]
        T % Temperature
        E % Potential at the end of collector
        
        sigmaeff % effective solid conductivity

        % Mesh properties (setup before simulation) s
        X
        Xb
        N
        dombin
                
    end
    
    methods
        function obj = currentCollector(compname, T, dims, grid)
            
            obj = obj@FvModel(compname);
            
            obj.am = currentCollectorAM(T);
            obj.am.eps = 1;
            
            obj.T = T;

            obj.Grid = genSubGrid(G, cells);
            
            obj.varnames = {'phi'};
            obj.varsizes = obj.Grid.cells.num;
        end
        
    end
end

