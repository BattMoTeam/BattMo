function [] = animateVoltage(model, t, states_ccpe, filename)

background =  [0, 60, 101]./255;%0.2 .* [1 1 1]; %


%% Plot the electric potential in the electrodes
h = figure;

for i = 1:10:length(t)
    disp(t(i));
    plot(t(1:i)./3600, E(1:i), '-w', 'LineWidth', 5);
    
    xlabel('Time  /  h')
    ylabel('Voltage  /  V')
    
    
ax = gca;
ax.XColor = 'w';

ax.YColor = 'w';

x0=5;
y0=5;
width=20;
height=20;
set(gcf,'units','centimeter','position',[x0,y0,width,height], 'color',background);
set(gca, 'color', background, 'FontSize', 14);

xmin = min(t)./3600;
xmax = max(t)./3600;

ymin = min(E);
ymax = max(E);

ylim([ymin, ymax]);
xlim([xmin, xmax]);

    
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
