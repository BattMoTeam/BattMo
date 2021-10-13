function jsonfd = pickField(jsonstruct, fdname)
    if ~isempty(jsonstruct) 
        if isfield(jsonstruct, fdname)
            jsonfd = jsonstruct.(fdname);
        end
    else
        jsonfd = [];
    end
end

