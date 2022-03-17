function out = blueorange(M)
%BLUEORANGE A divergent colormap from blue to orange
%   out = blueorange()      returns the colormap with the default
%                           resolution.  
%   out = blueorange(M)     returns an M-by-3 matrix containing a
%                           colormap. The colormap linearly interpolated to
%                           M discrete levels.
%
%   raw data source: https://sciviscolor.org/home/colormaps/contrasting-divergent-colormaps/

raw = [ 0.0863    0.0039    0.2980;
        0.1137    0.0235    0.4510;
        0.1059    0.0510    0.5098;
        0.0392    0.0392    0.5608;
        0.0314    0.0980    0.6000;
        0.0431    0.1647    0.6392;
        0.0549    0.2431    0.6784;
        0.0549    0.3176    0.7098;
        0.0510    0.3961    0.7412;
        0.0392    0.4667    0.7686;
        0.0314    0.5373    0.7882;
        0.0314    0.6157    0.8118;
        0.0235    0.7098    0.8314;
        0.0510    0.8000    0.8510;
        0.0706    0.8549    0.8706;
        0.2627    0.9020    0.8627;
        0.4235    0.9412    0.8745;
        0.5725    0.9647    0.8353;
        0.6588    0.9804    0.8431;
        0.7647    0.9804    0.8667;
        0.8275    0.9804    0.8863;
        0.8902    0.9882    0.9255;
        0.9137    0.9882    0.9373;
        1.0000    1.0000    0.9725;
        0.9882    0.9882    0.9059;
        0.9922    0.9725    0.8039;
        0.9922    0.9647    0.7137;
        0.9882    0.9569    0.6431;
        0.9804    0.9176    0.5098;
        0.9686    0.8745    0.4078;
        0.9490    0.8235    0.3216;
        0.9294    0.7765    0.2784;
        0.9098    0.7176    0.2353;
        0.8902    0.6588    0.1961;
        0.8784    0.6196    0.1686;
        0.8706    0.5490    0.1569;
        0.8510    0.4745    0.1451;
        0.8314    0.4118    0.1333;
        0.8118    0.3451    0.1137;
        0.7882    0.2667    0.0941;
        0.7412    0.1843    0.0745;
        0.6902    0.1255    0.0627;
        0.6196    0.0627    0.0431;
        0.5490    0.0275    0.0706;
        0.4706    0.0157    0.0902;
        0.4000    0.0039    0.1020;
        0.1882         0    0.0706  ];
    
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
