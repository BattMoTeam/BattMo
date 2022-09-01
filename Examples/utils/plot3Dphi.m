G = model.G;

ne    = 'NegativeElectrode';
pe    = 'PositiveElectrode';
eac   = 'ActiveMaterial';
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
    plotCellData(model.(ne).(am).G, state.(ne).(am).phi, plotparams{:});
    colorbar
    view(viewparam);
    title('phi (negative elde)');
    
    subplot(3, 2, 2);
    cla
    plotCellData(model.(pe).(am).G, state.(pe).(am).phi, plotparams{:});
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
    plot((time(1 : ind)/hour), E(1 : ind), '*-');
    xlabel('hours');
    ylabel('E');
    axis([0, max(time)/hour, min(E), max(E)])
    
    if dovideo
        frame = getframe(gcf);
        writeVideo(video, frame);
    end
    
    pause(0.1);
    
end

if dovideo
    close(video);
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
