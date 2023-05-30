function [flatjson, printFunction] = flattenJson(jsonstruct)

    flatjson = flattenJson_({}, jsonstruct, []);
    flatjson = reshape(flatjson, 2, [])';

    if nargout > 1
        printFunction = @(cellarray) printFunction_(cellarray);
    end

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


function printFunction_(cellarray)

    nrow = size(cellarray, 1);
    ncol = size(cellarray, 2);

    assert(ncol == 2, 'implementation only for 2 columns');
    
    strs = cell(nrow, ncol);
    ls   = zeros(ncol, 1);
    
    for ijson = 1 : nrow
        for col = 1 : ncol
            str = formattedDisplayText(cellarray{ijson, col});
            str = strtrim(str);
            str = regexprep(str, '\n', ',');
            ls(col) = max(strlength(str), ls(col));
            strs{ijson, col} = str;
        end
    end

    formatstr = sprintf('%%-%ds, %%-%ds\\n', ls(1), ls(2));
    for ijson = 1 : nrow
        fprintf(formatstr, strs{ijson, 1}, strs{ijson, 2});
    end

        
end
