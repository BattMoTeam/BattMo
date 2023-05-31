classdef FlatJsonViewer

    properties

        flatjson
        columnnames
        
    end

    methods

        function fjv  = FlatJsonViewer(flatjson, columnnames)

            fjv.flatjson = flatjson;

            if nargin < 2
                ncol = size(flatjson, 2);
                % default values for the column
                columnnames{1} = 'parameter name';
                if ncol == 2
                    columnnames{2} = 'parameter value';
                else
                    columnnames{2} = 'first set';
                    columnnames{3} = 'second set';
                    if ncol == 4
                        columnnames{4} = 'comparison';
                    end
                end
            end
            
            fjv.columnnames = columnnames;
            
        end


        function print(fjv, filterdesc)

            if nargin > 1
                fjv = fjv.filter(filterdesc);
            end
            
            cell2table(fjv.flatjson, 'VariableNames', fjv.columnnames)
            
        end

        function sortedfjv = sort(fjv, orderdesc)

            columnnames = fjv.columnnames;
            flatjson    = fjv.flatjson;

            ncol = numel(columnnames);

            if ischar(orderdesc)
                orderdesc = {orderdesc};
            end
            
            for iorder = numel(orderdesc) : -1 : 1

                r = regexprep(orderdesc{iorder}, ' +', '.*');
                ind = regexp(columnnames, r);
                ind = cellfun(@(res) ~isempty(res), ind);
                ind = find(ind);

                assert(numel(ind) == 1, 'regexp given is return too many match for column name');
                       
                [~, ia] = sort(flatjson(:, ind));
                flatjson = flatjson(ia, :);
                       
            end

            fjv.flatjson = flatjson;
            sortedfjv = fjv;

        end

        function filteredfjv = filter(fjv, filterdesc)

            if iscell(filterdesc{1})
                for ifilter = 1 : numel(filterdesc)
                    fjv = fjv.filter(filterdesc{ifilter});
                end
                filteredfjv = fjv;
                return
            end

            columnname  = filterdesc{1};
            filterval   = filterdesc{2};
            columnnames = fjv.columnnames;
            flatjson    = fjv.flatjson;

            r = regexprep(columnname, ' +', '.*');
            ind = regexp(columnnames, r);
            ind = cellfun(@(res) ~isempty(res), ind);
            ind = find(ind);
            assert(numel(ind) == 1);

            rowvals = flatjson(:, ind);
            r = regexprep(filterval, ' +', '.*');
            ind = regexp(rowvals, r);
            ind = cellfun(@(res) ~isempty(res), ind);
            ind = find(ind);

            flatjson = flatjson(ind, :);

            fjv.flatjson = flatjson;
            filteredfjv = fjv;

        end
        
    end
    
end
