function xv = removeAD(x)
%Remove AD state and compact 1 by n cell arrays to matrices
%
% SYNOPSIS:
%    v = removeAD(V);
%
% DESCRIPTION:
%   Removes AD variables
%
% REQUIRED PARAMETERS:
%   v   - Value to be converted.
%
% RETURNS:
%   v    - Value with no AD status.
%
% SEE ALSO:
%   double2ADI, double2GenericAD, ADI, GenericAD

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    if isnumeric(x) || islogical(x) || ischar(x)
        xv = x;
    elseif iscell(x)
        xv = applyFunction(@removeAD, x, 'uniformoutput', false);
    elseif isstruct(x)
        fn = fieldnames(x);
        for i = 1:numel(fn)
            f = fn{i};
            for j = 1:numel(x)
                x(j).(f) = removeAD(x(j).(f));
            end
        end
        xv = x;
    else
        try
            % If the object somehow has the property '.val' we use that. Note
            % that this is partially to fix an Octave bug when using function
            % handles on classes.
            xv = x.val;
        catch
            xv = x;
        end
    end
end
