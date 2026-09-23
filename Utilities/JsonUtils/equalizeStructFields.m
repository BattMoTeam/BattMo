function jsonstruct = equalizeStructFields(jsonstruct, fieldnamelists, varargin)
% same as equalizeStructField (without 's' at end) but iterate over the fieldnamelists;

    fielnamelist1 = fieldnamelists{1};

    for ifl = 2 : numel(fieldnamelists)
        
        fieldnamlist2 = fieldnamelists{ifl};

        jsonstruct = equalizeStructField(jsonstruct, fielnamelist1, fieldnamlist2, varargin{:});
        
    end

    
end

%{
Copyright 2021-2026 SINTEF Industry, Sustainable Energy Technology
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
