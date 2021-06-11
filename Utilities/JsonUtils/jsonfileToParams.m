function paramobj = jsonfileToParams(paramobj, filename)
    
    jsonsrc = fileread(filename);
    data = jsondecode(jsonsrc);

    paramobj = assignStructParams(paramobj, data);
    
end

