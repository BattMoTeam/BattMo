function loadModule(modulename)

    % Setup python environment
    if mrstPlatform('matlab')
        penv = pyenv;

        if penv.Version == ""

            % Try to setup python
            try

                pyenv('Version', '/usr/bin/python3.10');

            catch e

                warning('Cannot load module %s because no python installation is found. Try setting it manually, for example as pyenv(''Version'', ''/usr/bin/python3.10'');. You may have to install the correct libpython package (eg. libpython3.10) separately. To find the path to the executable you may use ''which python3''.\n', modulename);
                e.message
                e.ExceptionObject

            end

            % Try to add directory to the python path
            try

                dirname = fullfile(battmoDir(), 'Utilities', 'JsonUtils');

                if count(py.sys.path, dirname) == 0
                    dispif(mrstVerbose, 'Directory %s not found in py path, trying to insert it', dirname);

                    insert(py.sys.path, int32(0), dirname);
                    if count(py.sys.path, dirname) == 0
                        insert(py.sys.path, int64(0), dirname);
                    end
                    if count(py.sys.path, dirname) == 0
                        py.sys.path.insert(int32(0), dirname);
                    end
                    if count(py.sys.path, dirname) == 0
                        py.sys.path.insert(int64(0), dirname);
                    end

                    assert(count(py.sys.path, dirname))
                end

            catch e

                warning('Could not add directory to Python path. This may be due to an incompability between the MATLAB and Python versions, or that the correct libpython package is installed (eg. libpython3.10). See also https://se.mathworks.com/support/requirements/python-compatibility.html.');
                e.message
                e.ExceptionObject

            end
        end

    else
        warning('Calling python is not supported on this platform (%s). Only MATLAB is supported.', mrstPlatform('platform'));
    end

    % Try to load module
    try

        py.importlib.import_module(modulename);
        % mod = py.importlib.import_module(modulename);
        % py.importlib.reload(mod);

    catch e

        warning('Failed to load module %s', modulename);
        e.message
        e.ExceptionObject

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
