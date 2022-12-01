function fullfilename = resolveFileName(filename)
    
% If the filename is not an absolute path (starts with "/") or explicitly refer to the current directory (starts with
% "./"), then use battmo directory path.
    
    if strcmp(filename(1), '/') | strcmp(filename(1:2), './') |  strcmp(filename(1:3), '../')
        fullfilename = filename; 
    else
        fullfilename = fullfile(battmoDir(), filename);
    end
    
end

