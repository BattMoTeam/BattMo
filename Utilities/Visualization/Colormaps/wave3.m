function out = wave3(M)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here


raw = [ 0.843137, 0.894118, 0.960784;
        0.635294, 0.768627, 0.909804;
        0.482353, 0.678431, 0.858824;
        0.368627, 0.611765, 0.800000;
        0.278431, 0.556863, 0.729412;
        0.196078, 0.521569, 0.650980;
        0.160784, 0.513725, 0.600000;
        0.121569, 0.474510, 0.501961;
        0.098039, 0.439216, 0.439216;
        0.047059, 0.419608, 0.356863;
        0.094118, 0.431373, 0.290196;
        0.113725, 0.450980, 0.262745;
        0.145098, 0.478431, 0.254902;
        0.192157, 0.549020, 0.192157;
        0.309804, 0.619608, 0.247059;
        0.443137, 0.701961, 0.313725;
        0.584314, 0.780392, 0.388235;
        0.737255, 0.850980, 0.509804;
        0.878431, 0.878431, 0.615686;
        0.921569, 0.894118, 0.596078;
        0.901961, 0.835294, 0.521569;
        0.878431, 0.772549, 0.458824;
        0.870588, 0.705882, 0.384314;
        0.850980, 0.596078, 0.239216;
        0.800000, 0.529412, 0.254902;
        0.729412, 0.458824, 0.262745]; 
    
    if nargin == 0
        out = raw;
    else
        validateattributes(M,{'numeric'},{'nonnegative','scalar','real', 'integer'})
        s = size(raw);
        x = linspace(0,1,s(1));
        xo = linspace(0,1,M);
        out(:,1) = interp1(x,raw(:,1),xo);
        out(:,2) = interp1(x,raw(:,2),xo);
        out(:,3) = interp1(x,raw(:,3),xo);
    end

end




%{
Copyright 2021-2022 SINTEF Industry, Sustainable Energy Technology
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
