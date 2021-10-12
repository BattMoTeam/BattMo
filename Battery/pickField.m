function jsonfd = pickField(jsonstruct, fdname)
    if ~isempty(jsonstruct) 
        jsonstruct = resolveFileInputJson(jsonstruct);
        if isfield(jsonstruct, fdname)
            jsonfd = jsonstruct.(fdname);
            jsonfd = resolveFileInputJson(jsonfd);
        end
    else
        jsonfd = [];
    end
end

