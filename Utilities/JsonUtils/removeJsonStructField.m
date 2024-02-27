function jsonstruct = removeJsonStructField(jsonstruct, fdnames)

    fdname = fdnames{1};
    if numel(fdnames) == 1
        jsonstruct = rmfield(jsonstruct, fdname);
        return
    else
        jsonstruct.(fdname) = removeJsonStructField(jsonstruct.(fdname), fdnames(2:end));
    end
    
end

