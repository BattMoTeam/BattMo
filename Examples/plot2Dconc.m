G = model.G;

ne    = 'NegativeElectrode';
pe    = 'PositiveElectrode';
eac   = 'ElectrodeActiveComponent';
am    = 'ActiveMaterial';
cc    = 'CurrentCollector';
elyte = 'Electrolyte';
sep   = 'Separator';
thermal = 'ThermalModel';

h = figure(); 
set(h, 'Position', [107 28 1357 902]);

dovideo = false;

if dovideo
    filename = 'concentration.avi';
    video = VideoWriter(filename);
    video.FrameRate = 3;
    open(video);
end

for ind = 1 : numel(states)
    
    state = states{ind};
    
    figure(h);
    
    subplot(3, 2, 1);
    cla
    plotCellData(model.(elyte).G, state.(elyte).cs{1});
    colorbar
    title('cLi (elyte)');
    
    subplot(3, 2, 2);
    plot((time(1 : ind)/hour), Enew(1 : ind), '*-');
    xlabel('hours');
    ylabel('E');
    axis([0, max(time)/hour, min(Enew), max(Enew)])
    
    subplot(3, 2, 3);
    cla
    plotCellData(model.(ne).(eac).G, state.(ne).(eac).c);
    colorbar
    title('cLi (negative elde)');
    
    subplot(3, 2, 4);
    cla
    plotCellData(model.(pe).(eac).G, state.(pe).(eac).c);
    colorbar
    title('cLi (positive elde)');
    
    subplot(3, 2, 5);
    cla
    d = state.(ne).(eac).(am).cElectrode - state.(ne).(eac).c;
    plotCellData(model.(ne).(eac).G, d);
    colorbar
    title('negative elde : c(surface) - c(average)');
    
    subplot(3, 2, 6);
    cla
    d = state.(pe).(eac).(am).cElectrode - state.(pe).(eac).c;
    plotCellData(model.(pe).(eac).G, d);
    colorbar
    title('positive elde : c(surface) - c(average)');
    

    if dovideo
        frame = getframe(gcf);
        writeVideo(video, frame);
    end
    
    pause(0.1);
    
end

if dovideo
    close(video);
end

