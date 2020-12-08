classdef FvModel < handle
    
    properties
        Grid
        operators
        j_bcsource
    
        compname % name of the component
        
        varnames % names of the variables in the component
        varsizes % sizes of the variables in the component
        
    end
    
    function obj = FvModel(compname)
        obj.compname = compname;
    end
    
    function n = varnum(obj)
    % number of variables used for this model
        n = numel(varnames);
    end
    
    function n = N(obj)
    % number of grid cell (function used a short-cut)
        n = obj.Grid.cells.num;
    end
        
end

