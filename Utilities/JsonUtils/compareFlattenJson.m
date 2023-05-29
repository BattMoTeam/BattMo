function [flatjsondifferent, flatjsoncommon, flatjson1missing, flatjson2missing] = compareFlattenJson(flatjson1, flatjson2)

    fd1 = flatjson1(:, 1);
    fd2 = flatjson2(:, 1);

    [fd, ia, ib] = intersect(fd1, fd2);

    flatjsondifferent = {};
    flatjsoncommon    = {};
    
    assert(numel(ia) == numel(ib), 'Looks like I misunderstood something with intersect');
    
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

    flatjsoncommon    = vertcat(flatjsoncommon{:});
    flatjsondifferent = vertcat(flatjsondifferent{:});
    
    flatjson1missing = [];
    flatjson2missing = [];
    
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
