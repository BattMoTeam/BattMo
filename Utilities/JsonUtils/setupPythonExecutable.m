function ok = setupPythonExecutable(exec)

% Low level tool for setting the Python executable

    if pyenv().Version == ""

        if nargin == 0 || isempty(exec)

            % FIXME Base selection on MATLAB version and operating system.
            % https://se.mathworks.com/support/requirements/python-compatibility.html

            if isunix || ismac
                executables = {'/usr/bin/python3.11', ...
                               '/usr/bin/python3.10', ...
                               '/usr/bin/python3.9', ...
                               '/usr/bin/python3.8', ...
                               '/usr/bin/python3.7'};
            else % windows
                executables = {'C:\Python311', ...
                               'C:\Python310', ...
                               'C:\Python39', ...
                               'C:\Python38', ...
                               'C:\Python37'};
            end

            for k = 1:numel(executables)
                ok = setExec(executables{k});
                if ok
                    return;
                end
            end
        else
            ok = setExec(exec);
        end

    end

    assert(pyenv().Version ~= "", 'No valid Python executable set.');
    dispif(mrstVerbose, 'Python executable is %s\n', pyenv().Executable);

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
