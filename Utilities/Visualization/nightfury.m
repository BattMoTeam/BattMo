function nightfury(varargin)
% NIGHTFURY Sets the style of the current figure to the nightfury settings

aspect = 1./sqrt(2);
width = 29.7;
height = width * aspect;
discrete = [172, 223, 245;
    42, 174, 215;
    14, 120, 189;
    22, 59, 125;
    10, 9, 51;
    80, 14, 51;
    154, 54, 90;
    180, 62, 83;
    234, 101, 81;
    241, 133, 93;
    255, 171, 118;
    255, 199, 129;
    255, 232, 168] ./ 255;

map = lithiumIon.magma();

background = 0.2 .* [1 1 1];   

if strcmp(varargin, 'subplot')
                
    set(gca,'FontSize',14, ...
    'color',background, ...
    'ColorOrder', discrete)
    ax = gca;
    ax.XColor = 'w';
    ax.YColor = 'w';
    set(gcf,'units','centimeter',...
        'position',[1.5,1.5,width*1.6,height*1.3], ...
        'color',background)
    set(gca, 'FontName', 'Helvetica')

else
    set(gca,'FontSize',28, ...
        'color',background, ...
        'ColorOrder', discrete)
    ax = gca;
    ax.XColor = 'w';
    ax.YColor = 'w';
    set(gcf,'units','centimeter',...
        'position',[5,5,width,height], ...
        'color',background)
    set(gca, 'FontName', 'Helvetica')
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
