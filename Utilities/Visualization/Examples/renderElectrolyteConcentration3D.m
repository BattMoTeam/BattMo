close all
clc

filename = 'test.gif';

background = 0.2 .* [1 1 1]; % [0, 60, 101]./255;%

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
caxis([min(states{numStates}.elyte.cs{1}), max(states{numStates}.elyte.cs{1})]);

for ii = 1:numStates

p = plotCellData(model.elyte.G, states{ii}.elyte.cs{1});
p.LineStyle = 'none';

%light               % add a light
%lighting gouraud    % preferred lighting for a curved surface
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