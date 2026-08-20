function fn = setupParameters(fn, parameters)

    if isstruct(parameters) && numel(parameters) > 1
        parameters = num2cell(parameters);
    end
    
    for iparam = 1 : numel(parameters)

        param = parameters{iparam};
        
        fn.(param.name) = param.value;
        
    end

end
