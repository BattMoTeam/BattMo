thermal = 'ThermalModel';
G = model.G;


h = figure(); 
set(h, 'Position', [10 10 1700 500]);

dovideo = false;

if dovideo
    filename = 'temperature.avi';
    video = VideoWriter(filename);
    video.FrameRate = 3;
    open(video);
end

for ind = 1 : numel(states)

    subplot(1, 2, 1);
    plotCellData(G, states{ind}.(thermal).T)
    colorbar
    view([30, 32]);
    subplot(1, 2, 2);
    plot((time(1 : ind)/hour), Enew(1 : ind), '*-');
    xlabel('hours');
    ylabel('E');
    axis([0, max(time)/hour, min(Enew), max(Enew)])
    
    if dovideo
        frame = getframe(gcf);
        writeVideo(video, frame);
    end
    
    pause(0.1);
    
end

if dovideo
    close(video);
end

