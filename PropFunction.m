classdef PropFunction
    
    properties
        
        name
        fn
        namespace

    end
    
    methods
        
        function propfunction = PropFunction(name, fn, namespace)
            
            propfunction.name = name;
            propfunction.fn = fn;
            propfunction.namespace = namespace;
            
        end
        
    end
end
