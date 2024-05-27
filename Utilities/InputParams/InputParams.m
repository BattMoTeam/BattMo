classdef InputParams

    properties

        jsonstruct
        
    end
    
    methods

        function inputparams = InputParams(jsonstruct)

            inputparams = assignJsonParams(inputparams, jsonstruct);
            inputparams.jsonstruct = jsonstruct;
            
        end

        function inputparams = setParam(inputparams, names, val)
        % Note : same syntax as setProp in BaseModel
            if iscell(names) && (numel(names) > 1)
                name = names{1};
                names = names(2 : end);
                inputparams.(name) = setParam(inputparams.(name), names, val);
            elseif iscell(names) && (numel(names) == 1)
                name = names{1};
                if isnumeric(name)
                    inputparams{name} = val;
                else
                    inputparams.(name) = val;
                end
            else
                error('format not recognized');
            end
        end

        function var = getParam(inputparams, names)
        % Note : same syntax as getProp in BaseModel
            if iscell(names) && (numel(names) > 1)
                name = names{1};
                names = names(2 : end);
                var = getParam(inputparams.(name), names);
            elseif iscell(names) && (numel(names) == 1)
                name = names{1};
                if isnumeric(name)
                    var = inputparams{name};
                else
                    var = inputparams.(name);
                end
            else
                error('format not recognized');
            end
        end

        function jsonstruct = buildJsonStruct(inputparams);

            jsonstruct = inputparams.jsonstruct;
            
            inputparamsFds = propertynames(inputparams);

            for ind = 1 : numel(inputparamsFds)

                fd = inputparamsFds{ind};

                if isa(inputparams.(fd), 'InputParams')
                    subjsonstruct_fd = inputparams.(fd).buildJsonStruct();
                    clear subjsonstruct;
                    subjsonstruct.(fd) = subjsonstruct_fd;
                    jsonstruct = mergeJsonStructs({subjsonstruct, jsonstruct}, 'warn', false);
                end

            end

        end
        
        function inputparams = validateInputParams(inputparams)

        % Default automatic behaviour is that all the properties of inputparams that belong to class InputParams get
        % validated.

            inputparamsFds = propertynames(inputparams);

            for ind = 1 : numel(inputparamsFds)

                fd = inputparamsFds{ind};

                if isa(inputparams.(fd), 'InputParams')
                    inputparams.(fd) = inputparams.(fd).validateInputParams();
                end

            end

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
