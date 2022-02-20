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
        found = false;
        elt = list{ilist};
        imergelist = 1;
        while ~found & (imergelist <= numel(mergedlist))
            melt = mergedlist{imergelist};
            if eq(elt.varname, melt.varname)
                assert(melt.varname.dim == 1, 'cell are not yet supported');
                found = true;
                mergedlist{imergelist} = elt;
            end
            imergelist = imergelist + 1;            
        end
        if ~found
            mergedlist{end + 1} = elt;
        end

    end
    
end