function fn = setupParameters(fn, parameters)

    if isempty(parameters)
        return
    end
    
    if isstruct(parameters)
        if numel(parameters) > 1
            parameters = num2cell(parameters);
        else
            parameters = {parameters};
            fn = setupParameters(fn, parameters);
        end
    end
    
    for iparam = 1 : numel(parameters)

        param = parameters{iparam};
        
        fn.(param.name) = param.value;
        
    end

end
