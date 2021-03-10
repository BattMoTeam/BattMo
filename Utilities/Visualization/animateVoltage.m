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
