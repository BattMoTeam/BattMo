function mergedlist = mergeList(givenlist, list)
    
    if isempty(list)
        mergedlist = givenlist;
        return
    end
    
    if isempty(givenlist)
        mergedlist = list;
        return
    end
    
    % we only test type from first element
    testelt = list{1};

    if isa(testelt, 'VarName')
         mergedlist = mergeListVarNames(givenlist, list);
    elseif isa(testelt, 'PropFunction')
         mergedlist = mergeListPropFuncs(givenlist, list);
    else
        error('type not recognized');
    end

end

function mergedlist = mergeListVarNames(givenlist, list)
% We keep varname in place
    
    mergedlist = givenlist;
    
    for ilist  = 1 : numel(list)
        found = false;
        imergelist = 1;
        while ~found & imergelist <= numel(mergedlist)
            if eq(list{ilist}, mergedlist{imergelist})
                found = true;
            end
            imergelist = imergelist + 1;            
        end
        if ~found
            mergedlist{end + 1} = list{ilist};
        end
    end
    
end

function mergedlist = mergeListPropFuncs(givenlist, list)
% We overide the property in givenlist by a matching entry in list

    mergedlist = givenlist;
    
    for ilist  = 1 : numel(list)

        allfound = false;
        elt = list{ilist};
        imergedlist = 1;
        
        while ~allfound & (imergedlist <= numel(mergedlist))

            melt = mergedlist{imergedlist};
            res = decomposeVarName(elt.varname, melt.varname);
            
            if res.found
                
                newelt = elt;
                newelt.varname = res.intersecVarname;

                mergedlist = horzcat(mergedlist(1 : imergedlist - 1), {newelt}, mergedlist(imergedlist + 1 : end));
                
                if ~isempty(res.reminderGivenVarname)
                   newelt = melt;
                   newelt.varname = res.reminderGivenVarname;
                   imergedlist = imergedlist + 1;
                   mergedlist = horzcat(mergedlist(1 : imergedlist - 1), {newelt}, mergedlist(imergedlist + 1 : end));
                end
                
                if isempty(res.reminderVarname)
                    allfound = true;
                else
                    elt.varname = res.reminderVarname;
                end
                
            end
            
            imergedlist = imergedlist + 1;            
        end

        if ~allfound
            mergedlist{end + 1} = elt;
        end

    end
    
end