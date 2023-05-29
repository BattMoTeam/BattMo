function flatjson = flattenJson(jsonstruct)

    flatjson = flattenJson_({}, jsonstruct, []);
    flatjson = reshape(flatjson, 2, [])';

end

function flatjson = flattenJson_(flatjson, jsonstruct, prefix)

    dostruct = false;
    
    if isobject(jsonstruct)
        fds = properties(jsonstruct);
        dostruct = true;
    elseif isstruct(jsonstruct)    
        fds = fieldnames(jsonstruct);
        dostruct = true;
    else
        flatjson{end + 1} = prefix;
        flatjson{end + 1} = jsonstruct;
    end

    if dostruct
        for ifd = 1 : numel(fds)
            fd = fds{ifd};
            if isempty(prefix)
                subprefix = fd;
            else
                subprefix = sprintf('%s.%s', prefix, fd);
            end
            subjsonstruct = jsonstruct.(fd);
            if iscell(subjsonstruct) && numel(subjsonstruct) > 1
                subjsonstruct = {subjsonstruct};
            end
            flatjson = flattenJson_(flatjson, subjsonstruct, subprefix);
        end
    end
    
end

