function paramobj = jsonfileToParams(paramobj, filename)

    fileroot = '../..';
    fullfilename = fullfile(fileroot, filename);
    jsonsrc = fileread(filename);
    data = jsondecode(jsonsrc);

    paramobj = assignStructParams(paramobj, data);
    
end

