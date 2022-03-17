ne      = 'NegativeElectrode';
pe      = 'PositiveElectrode';
eac     = 'ActiveMaterial';
am      = 'ActiveMaterial';
thermal = 'ThermalModel';
elyte   = 'Electrolyte';
cc      = 'CurrentCollector';

Ts = cellfun(@(state) state.(thermal).T, states, 'uniformoutput', false);
G = model.G;

if strcmp(modelcase, '1D')

    thermal = 'ThermalModel';

    tfields = {'T', ...
               'jHeatOhmSource', ...
               'jHeatBcSource', ...
               'jHeatChemicalSource', ...
               'jHeatReactionSource'};
    
    mnames = {{'Electrolyte'}, ...
              {'PositiveElectrode','ActiveMaterial'}, ...
              {'NegativeElectrode','ActiveMaterial'}, ...
              {'NegativeElectrode','CurrentCollector'}, ...
              {'PositiveElectrode','CurrentCollector'}};    
    
    doFixedTempScale = false;
    if doFixedTempScale
        tM = max(states{1}.(thermal).T);
        tm = min(states{1}.(thermal).T);
        for i = 1 : numel(states)
            tM = max(tM, max(states{i}.(thermal).T));
            tm = min(tm, min(states{i}.(thermal).T));
        end
        tM = tM + 1e-1*(tM - tm);
        tm = tm - 1e-1*(tM - tm);
    end
    
    for i = 1 : numel(states)
        
        state = states{i}; 

        for k = 1 : numel(ffields)

            ffield = ffields{k};
            
            figure(h1)
            subplot(2, 3, k)
            cla, hold on
            
            if (strcmp(ffield, 'phi') || strcmp(ffield, 'j'))
                inds = 1 : 5; 
            else
                inds = 1 : 3; 
            end
            
            for ind = inds
                
                mname = mnames{ind}; 
                submodel = model.getSubmodel(mname); 
                substate = model.getProp(state, mname); 
                
                if strcmp(ffield, 'c') && (ind == 1)
                    var = substate.cs{1}; 
                else
                    var = substate.(ffield); 
                end
                
                if (k < 3)
                    plot(submodel.G.cells.centroids(:, 1), var, '* - ')
                else
                    iface = all(submodel.G.faces.neighbors > 0, 2); 
                    plot(submodel.G.faces.centroids(iface, 1), var, '* - ')
                end
                
                subtitle(ffield)
                
            end
        end
     
        % plot temperature
        subplot(2, 3, 5);
        cc = model.G.cells.centroids(:, 1);
        vols = model.G.cells.volumes;
        plot(cc, state.(thermal).T, '* - ');
        if doFixedTempScale
            axis([min(c), max(c), tm, tM]);
        end
        
        for k = 1 : numel(tfields)
            
            tfield = tfields{k};

            figure(h2)
            subplot(2, 3, k)
            cla, hold on
            val = state.(thermal).(tfield);
            if contains(tfield, 'Source')
                % we divide with volume if we have a source term.
                val = val./vols;
            end
            plot(cc, val, '* - ');
            subtitle(tfield);
            
        end
        
        drawnow;
        pause(0.01);
            
    end
    
end 


close all

dovideo = false;

if dovideo
    filename = 'temperature.avi';
    video = VideoWriter(filename);
    video.FrameRate = 3;
    open(video);
end

tM = max(Ts{1});
tm = min(Ts{1});
for i = 1 : numel(Ts)
    tM = max(tM, max(Ts{i}));
    tm = min(tm, min(Ts{i}));
end
tM = tM + 1e-1*(tM - tm);
tm = tm - 1e-1*(tM - tm);
        
for i = 1 : numel(Ts)
    clf        
    plotCellData(G, Ts{i}, 'edgealpha', 0.1)
    caxis([tm, tM]);
    titlestr = sprintf('Temperature, time = %g hour', states{i}.time/hour);
    title(titlestr)
    colorbar
    if dovideo
        frame = getframe(gcf);
        writeVideo(video, frame);
    end
    drawnow
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
