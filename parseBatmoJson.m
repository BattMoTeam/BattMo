function jsonstruct = parseBatmoJson(filename)
    filename = fullfile(batmoDir(), filename);
    jsonsrc = fileread(filename);
    jsonstruct = jsondecode(jsonsrc);
    jsonstruct = resolveFileInputJson(jsonstruct);
end
