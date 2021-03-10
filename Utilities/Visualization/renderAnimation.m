fps = 20;
background = 0.2 .* [1 1 1]; % [0, 60, 101]./255;%


%% Plot the electric potential in the electrodes
h = figure;
filename = 'neElectricPotential.gif';
colormap(blueorange(1000));
ax = gca;
ax.XColor = 'w';
ax.YColor = 'w';
x0=5;
y0=5;
width=20;
height=20;
set(gcf,'units','centimeter','position',[x0,y0,width,height], 'color',background);
set(gca, 'color',background, 'FontSize', 14);
caxis([0, 4.2]);
c = colorbar;
c.Color = 'w';

Gccne = model.ccne.G;
Gne = model.ne.G;

Gccpe = model.ccpe.G;
Gpe = model.pe.G;

xmin = min(Gccne.nodes.coords(:,1));
xmax = max(Gccpe.nodes.coords(:,1));

ymin = min(Gccne.nodes.coords(:,2));
ymax = max(Gne.nodes.coords(:,2));

ylim([ymin, ymax]);
xlim([xmin, xmax]);

for i = 1:10:length(states_ccpe)
    disp(t(i));
    p1 = plotCellData(Gccne, states_ccne{i}.phi);
    p1.EdgeColor = 'none';
    hold on
    p2 = plotCellData(Gne, states_ne{i}.phi);
    p2.EdgeColor = 'none';
    
    p3 = plotCellData(Gccpe, states_ccpe{i}.phi);
    p3.EdgeColor = 'none';
    p4 = plotCellData(Gpe, states_pe{i}.phi);
    p4.EdgeColor = 'none';
        drawnow
   
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
