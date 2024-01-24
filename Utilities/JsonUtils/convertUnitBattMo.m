function [val, isConverted] = convertUnitBattMo(val)
% We check if the field is a numerical parameter entry with units, that
% val has two fields
% - value : The numerical value
% - unit : A string with unit that can be evaluated in BattMo.
%          we use MRST support for unit, see "battmoDir()/Externals/mrst/mrst-core/utils/units/"
%          An example is "ampere/(centi*meter)^2"

    isConverted = false;

    if isfield(val, 'value') && isfield(val, 'unit')
        if numel(val) == 1
            isConverted = true;
            % This is a numerical parameter that requires a unit conversion
            if ~isempty(val.unit)
                str = sprintf('val = %g*%s;', val.value, val.unit);
                eval(str);
            else
                % Empty unit. Assume value is given in SI units
                val = val.value;
            end
        else
            vals = [];
            for ival = 1 : numel(val)
                cval = convertUnitBattMo(val(ival));
                vals(end + 1) = cval;
            end
            isConverted = true;
            val = vals;
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
