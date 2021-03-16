close all
clc

filename = 'renderPotential3D.gif';


background = [0, 60, 101]./255;%

h = figure;
ax = gca;
ax.XColor = 'w';
ax.YColor = 'w';
x0=5;
y0=5;
width=27;
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
caxis([0, 4.2]);

for ii = 1:numStates

pltelyte = plotCellData(model.elyte.G, states{ii}.elyte.phi);
hold on
pltccpe = plotCellData(model.ccpe.G, states{ii}.ccpe.phi);
pltpe = plotCellData(model.pe.G, states{ii}.pe.am.phi);
pltne = plotCellData(model.ne.G, states{ii}.ne.am.phi);
pltccne = plotCellData(model.ccne.G, states{ii}.ccne.phi);
pltelyte.LineStyle = 'none';
pltccpe.LineStyle = 'none';
pltpe.LineStyle = 'none';
pltne.LineStyle = 'none';
pltccne.LineStyle = 'none';

%light               % add a light
%lighting gouraud    % preferred lighting for a curved surface
material dull
view(0,-90)
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