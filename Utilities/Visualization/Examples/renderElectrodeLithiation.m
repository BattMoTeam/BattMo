close all
clc

filename = 'renderElectrodeLithiation.gif';

background = [0, 60, 101]./255;%

h = figure;
ax = gca;
ax.XColor = 'w';
ax.YColor = 'w';
x0=5;
y0=5;
width=1.5*27;
height=21.9;
set(gcf,'units','centimeter','position',[x0,y0,width,height], 'color',background);
set(gca, 'color',background, 'FontSize', 14);

c = colorbar;
c.Color = 'w';
cmap = flipud(orangeblue(1000));
colormap(cmap);

[maxDims] = max(model.G.nodes.coords);
d0 = max(maxDims);
pbaspect([maxDims ./ d0]);
numStates = length(find(~cellfun('isempty', states)));
caxis([0,1]);

for ii = 1:numStates

subplot(1,2,1)
pbaspect([maxDims ./ d0]);
p = plotCellData(model.pe.G, states{ii}.pe.am.Li ./ model.ne.am.Li.cmax);
p.LineStyle = 'none';
title('Positive Electrode', 'Color', 'w', 'FontSize', 24)

%light               % add a light
%lighting gouraud    % preferred lighting for a curved surface
if ii == 1
    camzoom(2)
end
material dull
view(60,30)
axis off      % remove axis

subplot(1,2,2)
pbaspect([maxDims ./ d0]);
p = plotCellData(model.ne.G, states{ii}.ne.am.Li ./ model.ne.am.Li.cmax);
p.LineStyle = 'none';
title('Negative Electrode', 'Color', 'w', 'FontSize', 24)

%light               % add a light
%lighting gouraud    % preferred lighting for a curved surface
if ii == 1
    camzoom(2)
end
material dull
view(60,30)
axis off      % remove axis

drawnow
frame = getframe(h); 
im = frame2im(frame); 
[imind,cm] = rgb2ind(im,256); 
if ii == 1
    imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
else 
   imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0); 
end
pause(0.1)



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
