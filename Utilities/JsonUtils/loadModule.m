function loadModule(modulename)

    if nargin == 0
        modulename = 'validationJsonScript';
    end

    % Setup python
    if mrstPlatform('matlab')
        pe = pyenv;

        if pe.Version == ""

            % Try to setup python
            try
                pyenv('Version', '/usr/bin/python3.10');
            catch
                warning('Cannot load module %s because no python installation is found. Try setting it manually, for example as pyenv(''Version'', ''/usr/bin/python3.10'');. You may have to install the correct libpython package (eg. libpython3.10) separately. To find the path to the executable you may use ''which python3''.\n', modulename);
            end

            % Try to add the path
            try
                rootdirname = fileparts(mfilename('fullpath'));
                dirname = fullfile(rootdirname);
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

    % Try to load module
    try
        mod = py.importlib.import_module(modulename);
        py.importlib.reload(mod);
    catch
        warning('Failed to load module %s', modulename);
    end

end



%{
Copyright 2021-2023 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
