function setupPythonPath(dirs)

% Low level tool for adding directories to the Python path. The
% directories should be full paths.

    if nargin == 0
        dirs = fullfile(battmoDir(), 'Utilities', 'JsonUtils');
    end

    if ~iscell(dirs)
        dirs = {dirs};
    end

    for k = 1:numel(dirs)

        dirname = dirs{k};

        if count(py.sys.path, dirname) == 0

            dispif(mrstVerbose, 'Directory %s is not found in the Python path, trying to adding it...\n', dirname);

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

        end

        assert(count(py.sys.path, dirname), 'Unable to add directory %s to the Python path');

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
