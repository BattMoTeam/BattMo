classdef couplingTerm
    
    properties
        name
        componentnames
        couplingcells
        couplingfaces
    end
    
    methods
        
        function obj = couplingTerm(name, compnames)
            obj.name = name;
            obj.compnames = compnames;
        end
        
    end
    
end
