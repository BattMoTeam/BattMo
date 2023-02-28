function loadModule(modulename)

    if nargin == 0
        modulename = 'validationJsonScript';
    end

    % Setup python
    if mrstPlatform('matlab')
        pe = pyenv;
        if pe.Version == ""
            warning('Cannot load module %s because no python installation is found. Try setting it manually, for example as pyenv(''Version'', ''/usr/bin/python3.10'');. You may have to install the correct libpython package (eg. libpython3.10) separately. To find the path to the executable you may use ''which python3''.\n', modulename);
        else
            try
                rootdirname = fileparts(mfilename('fullpath'));
                dirname = fullfile(rootdirname, '..', '..', 'Utilities', 'JsonUtils');
                pypath = cell(py.sys.path);
                found = false;
                for k = 1:numel(pypath)
                    if strcmpi(char(pypath{k}), dirname)
                        found = true;
                    end
                end
                if ~found
                    insert(py.sys.path, int32(0), dirname);
                end
            catch
                warning('Could not add directory to Python path. This may be due to an incompability between the MATLAB and Python versions, or that the correct libpython package is installed (eg. libpython3.10). See also https://se.mathworks.com/support/requirements/python-compatibility.html.');
            end
        end
    else
        warning('Calling python is not supported on this platform (%s). Only MATLAB is supported.', mrstPlatform('platform'));
    end

    mod = py.importlib.import_module(modulename);
    py.importlib.reload(mod);

end
