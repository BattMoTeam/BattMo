G = model.G;

ne    = 'NegativeElectrode';
pe    = 'PositiveElectrode';
eac   = 'ElectrodeActiveComponent';
cc    = 'CurrentCollector';
elyte = 'Electrolyte';
sep   = 'Separator';
thermal = 'ThermalModel';

h = figure(); 
set(h, 'Position', [10 10 900 1000]);

dovideo = false;

if dovideo
    filename = 'phi.avi';
    video = VideoWriter(filename);
    video.FrameRate = 3;
    open(video);
end

viewparam = [30, 47];
plotparams = {'linewidth', 0.01};

for ind = 1 : numel(states)
    
    state = states{ind};
    
    figure(h);
    
    subplot(3, 2, 1);
    cla
    plotCellData(model.(ne).(eac).G, state.(ne).(eac).phi, plotparams{:});
    colorbar
    view(viewparam);
    title('phi (negative elde)');
    
    subplot(3, 2, 2);
    cla
    plotCellData(model.(pe).(eac).G, state.(pe).(eac).phi, plotparams{:});
    colorbar
    view(viewparam);
    title('phi (positive elde)');
    
    subplot(3, 2, 3);
    cla
    plotCellData(model.(ne).(cc).G, state.(ne).(cc).phi, plotparams{:});
    colorbar
    view(viewparam);
    title('phi (negative cc)');
    
    subplot(3, 2, 4);
    cla
    plotCellData(model.(pe).(cc).G, state.(pe).(cc).phi, plotparams{:});
    colorbar
    view(viewparam);
    title('phi (positive cc)');
    
    subplot(3, 2, 5);
    cla
    plotCellData(model.(elyte).G, state.(elyte).phi, plotparams{:});
    colorbar
    view(viewparam);
    title('phi (elyte)');
    
    subplot(3, 2, 6);
    cla
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

