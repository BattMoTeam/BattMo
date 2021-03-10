function [] = animateLithiation(model, t, states_ne, states_pe, filename)

background =  [0, 60, 101]./255;% 0.2 .* [1 1 1]; %

%% Plot the electric potential in the electrodes
h = figure;
colormap(blueorange(1000));
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
caxis([0, 1]);
c = colorbar;
c.Color = 'w';

Gccne = model.ccne.G;
Gne = model.ne.G;

Gccpe = model.ccpe.G;
Gpe = model.pe.G;

xmin = min(Gne.nodes.coords(:,1));
xmax = max(Gpe.nodes.coords(:,1));

ymin = min(Gne.nodes.coords(:,2));
ymax = max(Gne.nodes.coords(:,2));

ylim([ymin, ymax]);
xlim([xmin, xmax]);

for i = 1:10:length(t)
    p1 = plotCellData(Gne, states_ne{i}.Li ./ model.ne.am.Li.cmax);
    p1.EdgeColor = 'none';
    hold on

    p2 = plotCellData(Gpe, states_pe{i}.Li ./ model.pe.am.Li.cmax);
    p2.EdgeColor = 'none';
        drawnow
    if nargin >= 5 
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
