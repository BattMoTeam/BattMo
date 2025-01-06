function loadModule(modulenames, varargin)

% Top level script for loading Python modules and (optionally) setup
% custom Python executable and custom Python paths).

    opt = struct('setupPython', true, ...
                 'dir', fullfile(battmoDir(), 'Utilities', 'JsonUtils'), ...
                 'exec', '', ...
                 'reload', false);
    opt = merge_options(opt, varargin{:});

    if mrstPlatform('matlab')

        if opt.setupPython
            setupPythonExecutable(opt.exec);
            setupPythonPath(opt.dir);
        end

        if ~iscell(modulenames)
            modulenames = {modulenames};
        end

        for k = 1:numel(modulenames)
            modulename = modulenames{k};
            try
                dispif(mrstVerbose, 'Loading module %s\n', modulename);
                py.importlib.import_module(modulename);

                if opt.reload

                    sys = py.importlib.import_module('sys');
                    if isfield(sys.modules, modulename)
                        remove(sys.modules, modulename);
                    else
                        warning('Module %s was not loaded before, so it cannot be reloaded.', modulename);
                    end

                    py.importlib.import_module(modulename);

                end

            catch e
                disp(e);
                error('Failed to load module %s', modulename);
            end
        end

    else
        warning('Calling Python is not supported on this platform (%s). Only MATLAB is supported.', mrstPlatform('platform'));
    end

end


%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
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
