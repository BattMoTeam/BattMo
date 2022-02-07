function dir = batmoDir()
    dir = fileparts(mfilename('fullpath'));
    dir = fullfile(dir, '..');
end
