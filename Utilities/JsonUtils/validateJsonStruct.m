function is_valid = validateJsonStruct(jsonstruct, varargin)

    opt = struct('useTmpFile', true, ...
                 'schema', fullfile(battmoDir(), 'Utilities', 'JsonSchemas', 'Simulation.schema.json'));

    [opt, extra] = merge_options(opt, varargin{:});

    if opt.useTmpFile

        % Write the json struct to a temporary file
        tempfilename = writeJsonStruct(jsonstruct);

        % Validate
        is_valid = validateJsonFiles({tempfilename});

    else

        jsonstr = jsonencode(jsonstruct);

        loadModule({'json', 'jsonschema'}, extra{:});

        fid = fopen(opt.schema, 'r');
        schemastr = fread(fid, '*char')';
        fclose(fid);

        pyschema = py.json.loads(schemastr);
        pydict = py.json.loads(jsonstr);

        try
            py.jsonschema.validate(pydict, pyschema);
            is_valid = true;
        catch ME
            disp(ME.message);
            is_valid = false;
        end

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
