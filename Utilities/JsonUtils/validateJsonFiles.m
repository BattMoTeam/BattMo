function is_valid = validateJsonFiles(jsonfiles)

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
        dispif(mrstVerbose, 'Validating %s\n', jsonfile);
        is_valid{k} = py.(modname).validate(jsonfile); %#ok
        assert(is_valid{k}, 'jsonfile %s is not valid', jsonfile);
    end

    if numel(jsonfiles) == 1
        is_valid = is_valid{1};
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
