function isValid = validateJsonFiles(jsonfiles)

    if ~iscell(jsonfiles)
        jsonfiles = {jsonfiles};
    end

    % Load module
    % FIXME Don't load the module every time
    modname = 'validateJsonFiles';
    loadModule(modname);

    % Validate using python script
    for k = 1:numel(jsonfiles)
        jsonfile = jsonfiles{k};
        dispif(mrstVerbose, '\nValidating %s\n\n', jsonfile);
        try
            isValid{k} = py.(modname).validate(battmoDir(), jsonfile, pyargs('verbose', mrstVerbose)); %#ok
        catch e
            disp(e);
            error('Error when running the validation');
        end
        assert(isValid{k}, 'jsonfile %s is not valid', jsonfile);
    end

    if numel(jsonfiles) == 1
        isValid = isValid{1};
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
