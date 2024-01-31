function inputparams = assignStructParams(inputparams, structdata)
    
    fields_sd = fieldnames(structdata);
    fields_pobj = fieldnames(inputparams);
    
    for ind = 1 : numel(fields_sd)
        
        fd = fields_sd{ind};
        
        % if isclass(inputparams)
            % if inputparams is a class, we check here that it the field fd matches a property of the class
            % assert(ismember(fd, fields_pobj), 'field in input data is not recognized');
        % end
        
        if isstruct(structdata.(fd)) && isfield(structdata.(fd), 'isFile') && structdata.(fd).isFile
            filename = structdata.(fd).filename;
            inputparams.(fd) = jsonfileToParams(inputparams.(fd), filename);
        elseif isnumeric(structdata.(fd))
            inputparams.(fd) = structdata.(fd);
        elseif ischar(structdata.(fd))
            inputparams.(fd) = structdata.(fd);
        elseif iscell(structdata.(fd))
            inputparams.(fd) = structdata.(fd);
        elseif isstruct(structdata.(fd))
            inputparams.(fd) = assignStructParams(inputparams.(fd), structdata.(fd));
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
