function sum = sha1sum(varargin)

    if mrstPlatform('matlab')

        % Possibly make this persistent
        md = java.security.MessageDigest.getInstance('SHA-1');

        % Serialize the input arguments to byte arrays
        bytes = cellfun(@(x) getByteStreamFromArray(x), varargin, 'UniformOutput', false);

        % Sort to get consistent hash, i.e. sha1sum('a', 'b') ==
        % sha1sum('b', 'a')
        flat = sort([bytes{:}]);

        % Hash
        md.update(flat);
        hash = typecast(md.digest(), 'uint8');

        % Convert to lowercase hex string
        sum = lower(reshape(dec2hex(hash), 1, []));

    else

        error('Only supported on MATLAB platform');

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
