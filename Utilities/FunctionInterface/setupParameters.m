function fn = setupParameters(fn, parameters)

    if isempty(parameters)
        return
    end

    fdnames = fieldnames(parameters);
    for iparam = 1 : numel(fdnames)

        fdname = fdnames{iparam};
        fn.(fdname) = parameters.(fdname);
        
    end

end
