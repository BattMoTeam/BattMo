function I = rampupControl(t, tup, inputI, varargin)
%
%
% SYNOPSIS:
%   function I = rampupControl(t, tup, inputI, varargin)
%
% DESCRIPTION: Generate a function that gradually increase the output I to the given inputI over the period given by tup
%
% PARAMETERS:
%   t        - current time value
%   tup      - rampup time: The value is increased in the interval from zero to tup
%   inputI   - Target value after rampup
%   varargin - used to set the rampup profile (sineup is default, can be set to linear)
%
% RETURNS:
%   I - value returned for the current time (the other values being given)
%


    opt = struct('rampupcase', 'sineup');
    opt = merge_options(opt, varargin{:});
   
    rampupcase = opt.rampupcase;
    
    switch rampupcase
        
      case 'sineup'
        I = (t <= tup) .* sineup(0, inputI, 0, tup, t) + (t > tup) .* inputI;
        
      case 'linear'
        I = (t <= tup).*t./tup .* inputI +  (t > tup) .* inputI;
        
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
