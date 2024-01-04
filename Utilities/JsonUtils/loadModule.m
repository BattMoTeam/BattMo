function loadModule(modulename, varargin)

% Top level script for loading Python modules and (optionally) setup
% custom Python executable and custom Python paths).

    opt = struct('setupPython', true, ...
                 'dir', fullfile(battmoDir(), 'Utilities', 'JsonUtils'));
    [opt, extra] = merge_options(opt, varargin{:});

    if mrstPlatform('matlab')

        if opt.setupPython
            setupPythonExecutable(extra{:});
            setupPythonPath(opt.dir);
        end

        try
            py.importlib.import_module(modulename);
        catch e
            disp(e);
            error('Failed to load module %s', modulename);
        end

    else
        warning('Calling Python is not supported on this platform (%s). Only MATLAB is supported.', mrstPlatform('platform'));
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
