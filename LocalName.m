classdef LocalName < VarName
    
    properties
    end
    
    methods
        function varname = LocalName(name, index)
            varname = varname@VarName('', name, index);
        end
        
        function varname = AttachLocalName(varname, model)
            varname.namespace = model.namespace;
        end
    end
end