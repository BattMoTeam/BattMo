function [flatjsondifferent, flatjsoncommon, flatjsonmissing, printFunction] = compareFlattenJson(flatjson1, flatjson2)

    fd1 = flatjson1(:, 1);
    fd2 = flatjson2(:, 1);

    [fd, ia, ib] = intersect(fd1, fd2);

    flatjsondifferent = {};
    flatjsoncommon    = {};
    flatjsonmissing   = {};
    
    for ii = 1 : numel(ia)
        val1 = flatjson1{ia(ii), 2};
        val2 = flatjson2{ib(ii), 2};
        isequal = compareValue(val1, val2);
        entry = {fd{ii}, val1, val2};
        if isequal
            flatjsoncommon{end + 1} = entry;
        else
            flatjsondifferent{end + 1} = entry;
        end
        
    end

    ismissing1 = true(numel(fd2), 1);
    ismissing1(ib) = false;
    ismissing1 = find(ismissing1);
    for ii = 1 : numel(ismissing1)
        entry = {flatjson2{ismissing1(ii), 1}, NaN, flatjson2{ismissing1(ii), 2}};
        flatjsonmissing{end + 1} = entry;
    end
    
    ismissing2 = true(numel(fd1), 1);
    ismissing2(ia) = false;
    ismissing2 = find(ismissing2);
    for ii = 1 : numel(ismissing2)
        entry = {flatjson1{ismissing2(ii), 1}, flatjson1{ismissing2(ii), 2}, NaN};
        flatjsonmissing{end + 1} = entry;
    end

    flatjsoncommon    = vertcat(flatjsoncommon{:});
    flatjsondifferent = vertcat(flatjsondifferent{:});
    flatjsonmissing   = vertcat(flatjsonmissing{:});

    printFunction = @(jsondiff) printFunction_(jsondiff);
end


function isequal = compareValue(val1, val2)
    
    if isstruct(val1) && isstruct(val2)
        % val1 and val2 are structs
        fds1 = fieldnames(val1);
        fds2 = fieldnames(val2);
        fds = intersect(fds1, fds2);
        if numel(fds) == numel(fds1) && numel(fds) == numel(fds2)
            % val1 and val2 have the same fields
            for ifd = 1 : numel(fds)
                % We compare the values of each field
                fd = fds{ids}
                subisequal = compareValue(val1.(fd), val2.(fd));
                if ~subisequal
                    isequal = false;
                    return
                end
            end
            isequal = true;
        else
            isequal = false
        end
        return;
    end
    
    if iscell(val1) && iscell(val2)
        % val2 and val2 are cells
        if numel(val1) == numel(val2)
            % we compare the value of each cell
            for ival = 1 : numel(val1)
                subisequal = compareValue(val1{ival}, val2{ival});
                if ~subisequal
                    isequal = false;
                    return
                end
            end
            isequal = true;
        else
            isequal = false
        end
        return
    end
    
    try
        isequal = eq(val1, val2);
    catch
        isequal = false;
    end

end


function printFunction_(jsondiff)

    n = size(jsondiff, 1);
    strs = cell(n, 3);
    ls = zeros(3, 1);
    for ijson = 1 : n
        for col = 1 : 3
            str = formattedDisplayText(jsondiff{ijson, col});
            str = strtrim(str);
            str = regexprep(str, '\n', ',');
            ls(col) = max(strlength(str), ls(col));
            strs{ijson, col} = str;
        end
    end

    formatstr = sprintf('%%-%ds, %%-%ds, %%-%ds\\n', ls(1), ls(2), ls(3));

    for ijson = 1 : n
        fprintf(formatstr, strs{ijson, 1}, strs{ijson, 2}, strs{ijson, 3});
    end
end
