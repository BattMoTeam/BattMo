close all
clc

filename = 'stackCell.gif';
nGif = 20;

background = [0, 60, 101]./255;%
temp = copper(10);
colorNECC = temp(8,:);

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


[maxDims] = max(model.G.nodes.coords);
d0 = max(maxDims);
pbaspect([maxDims ./ d0]);

p = plotFaces(model.ccpe.G);
p.LineStyle = 'none';
p.FaceColor = 0.7 .* ones(1,3);

              % add a light
% lightangle(45,90)
% lighting gouraud    % preferred lighting for a curved surface
% material metal
view(60,30)
axis off      % remove axis

drawnow
for i = 1:nGif
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    if i == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else 
       imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0); 
    end
end
pause(1)

hold on

p = plotFaces(model.pe.G);
p.LineStyle = 'none';
p.FaceColor = 0.5 .* ones(1,3);
drawnow
for i = 1:nGif
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0); 
end
pause(1)

hold on

p = plotFaces(model.elyte.G);
p.LineStyle = 'none';
p.FaceColor = 'w';
drawnow
for i = 1:nGif
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0); 
end
pause(1)

hold on

p = plotFaces(model.ne.G);
p.LineStyle = 'none';
p.FaceColor = 0 .* ones(1,3);
drawnow
for i = 1:nGif
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0); 
end
pause(1)

hold on

p = plotFaces(model.ccne.G);
p.LineStyle = 'none';
p.FaceColor = colorNECC;
drawnow
for i = 1:nGif
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0); 
end
pause(1)

drawnow


