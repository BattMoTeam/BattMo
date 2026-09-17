function [I, alpha] = controlfunc(time, Imax, tswitch, T, varargin)

    opt = struct('order', 'alpha-first');
    opt = merge_options(opt, varargin{:});
    
    switch opt.order
      case 'alpha-first'
        if time <= tswitch
            I = 0;
            alpha = time/tswitch;
        else
            I = (time - tswitch)/(T - tswitch)*Imax;
            alpha = 1;
        end
      case 'I-first'
        if time <= tswitch
            I     = time/tswitch*Imax;
            alpha = 0;
        else
            I     = Imax;
            % alpha = 0;
            alpha = (time - tswitch)/(T - tswitch);
        end
      otherwise
        error('order not recognized');
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
