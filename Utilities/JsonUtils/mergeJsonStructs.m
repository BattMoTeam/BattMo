function jsonstruct = mergeJsonStructs(jsonstructs)
% - We call a json structure (abbreviated jsonstruct) a MATLAB structure of the form that are produced by the jsondecode command
% - The input jsonstructs is list jsonstruct
% - The command mergeJsonStructs merges recursively all the jsonstruct contained in the list jsonstructs
% - If two jsonstruct assign the same field, then an error is sent

    if numel(jsonstructs) == 1
        jsonstruct = jsonstructs{1};
        return;
    end
    
    jsonstruct1 = jsonstructs{1};
    jsonstruct2 = jsonstructs{2};
    if numel(jsonstructs) >= 2
        jsonstructrests = jsonstructs(3 : end);
    else
        jsonstructrests = {};
    end

    fds1 = fieldnames(jsonstruct1);
    fds2 = fieldnames(jsonstruct2);

    jsonstruct = jsonstruct1;

    for ifd2 = 1 : numel(fds2)
        fd2 = fds2{ifd2};
        
        if ~ismember(fd2, fds1)
            % ok, add the substructure
            jsonstruct.(fd2) = jsonstruct2.(fd2);
        else
            if isstruct(jsonstruct.(fd2)) && isstruct(jsonstruct2.(fd2))
                % we have to check the substructure
                subjsonstruct = mergeJsonStructs({jsonstruct.(fd2), jsonstruct2.(fd2)});
                jsonstruct.(fd2) = subjsonstruct;
            else
                error('parameters are assigned twice');
            end
        end
    end

    if ~isempty(jsonstructrests)
        jsonstruct = mergeJsonStructs({jsonstruct, jsonstructrests{:}});
    end
    
end
