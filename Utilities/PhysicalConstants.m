classdef PhysicalConstants
    %PHYSICALCONSTANTS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        R       % Ideal gas constant
        F       % Faraday constant
        c0      % Standard concentration
        molarVolumeLi   % Molar Volume of pure Li
    end
    
    methods
        function obj = PhysicalConstants()
            %PHYSICALCONSTANTS Construct an instance of this class
            %   Detailed explanation goes here
            obj.R = 8.31446261815324;
            obj.F = 96485.3329;
            obj.c0 = 1000;
            obj.molarVolumeLi = 9 * 1E-6;
        end
    end
end




%{
Copyright 2021-2023 SINTEF Industry, Sustainable Energy Technology
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
