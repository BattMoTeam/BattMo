function [val, isConverted] = convertUnitBattMo(val)
% We check if the field is a numerical parameter entry with units, that
% val has two fields
% - value : The numerical value
% - unit : A string with unit that can be evaluated in BattMo.
%          we use MRST support for unit, see "battmoDir()/mrst/mrst-core/utils/units/"
%          An example is "ampere/(centi*meter)^2"

    isConverted = false;
    
    if isfield(val, 'value') && isfield(val, 'unit')
        isConverted = true;
        % This is a numerical parameter that requires a unit conversion
        if ~isempty(val.unit)
            str = sprintf('val = %g*%s;', val.value, val.unit);
            eval(str);
        else
            % Empty unit. Assume value is given in SI units
            val = val.value;
        end
    end

end
