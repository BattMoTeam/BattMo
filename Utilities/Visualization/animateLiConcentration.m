function [] = animateLiConcentration(model, t, states_elyte, filename)

fps = 20;
background =  [0, 60, 101]./255;%0.2 .* [1 1 1]; %


%% Plot the electric potential in the electrodes
h = figure;
colormap(flipud(wave3(1000)));
ax = gca;
ax.XColor = 'w';
ax.XTick = [];
ax.YColor = 'w';
ax.YTick = [];
x0=5;
y0=5;
width=10;
height=20;
set(gcf,'units','centimeter','position',[x0,y0,width,height], 'color',background);
set(gca, 'color',background, 'FontSize', 14);
caxis([0.700, 1.200]);
c = colorbar;
c.Color = 'w';

Gelyte = model.elyte.G;

xmin = min(Gelyte.nodes.coords(:,1));
xmax = max(Gelyte.nodes.coords(:,1));

ymin = min(Gelyte.nodes.coords(:,2));
ymax = max(Gelyte.nodes.coords(:,2));

ylim([ymin, ymax]);
xlim([xmin, xmax]);

for i = 1:10:length(t)
    p1 = plotCellData(Gelyte, states_elyte{i}.Li ./ 1000);
    p1.EdgeColor = 'none';
    drawnow
   
    if nargin >= 4      
    % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
        if i == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
        else 
          imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0); 
        end 
    end
    
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
