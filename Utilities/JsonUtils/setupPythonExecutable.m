function ok = setupPythonExecutable(varargin)

% Low level tool for setting the Python executable.

    opt = struct('executable', '');
    opt = merge_options(opt, varargin{:});

    if isempty(opt.executable)

        % FIXME Base selection on MATLAB version and operating system.
        % https://se.mathworks.com/support/requirements/python-compatibility.html

        executables = {'/usr/bin/python3.10', '/usr/bin/python3.9'};

        for k = 1:numel(executables)
            ok = setExec(executables{k});
            if ok
                return;
            end
        end
    else
        ok = setExec(opt.executable);
    end

    assert(ok, 'Unable to set Python executable');

end

function ok = setExec(exec)

    dispif(mrstVerbose, 'Trying to set Python executable to %s...\n', exec);

    % Put in try/catch to allow for errors
    try
        pyenv('Version', exec);
        dispif(mrstVerbose, 'Setting Python executable to %s\n', pyenv().Executable);
        ok = true;
    catch
        dispif(mrstVerbose, 'Unable to set Python executable to %s\n', exec);
        ok = false;
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
