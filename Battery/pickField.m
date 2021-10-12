function fd = pickField(jsonstruct, fdname)
    if ~isempty(jsonstruct) 
        if isfield(jsonstruct, 'isFile')
            fileroot = batmoDir(); 
            filename = jsonstruct.filename;
            fullfilename = fullfile(fileroot, filename);
            jsonsrc = fileread(fullfilename);
            parsedjson = jsondecode(jsonsrc);
            jsonstruct = parsedjson; 
        end
        if isfield(jsonstruct, fdname)
            fd = jsonstruct.(fdname);
        end
    else
        fd = [];
    end
end

