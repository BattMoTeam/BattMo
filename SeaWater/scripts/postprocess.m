clear all

figure
hold on

directories =  {'seawater2', 'seawater3'};

c = parula(floor(5));
colororder(gca, c);
set(gca, 'linestyleorder', {'-', '-.'});

for  idir = 1 : numel(directories) 

    directory = directories{idir};
    [bp, simlist] = setupSimList(directory);

    mrstModule add ad-core

    newsimlist = simlist;

    newsimlist = bp.filterSimList(simlist, 'nucMax', 1, 'nparams', @(str) ismember(str, {'10;40;40', '10;20;40', '10;40;20', '80;80;80', '40;40;40', '30;30;30'}));

    bp.printSimList(newsimlist);

    

    for isim = 1 : numel(newsimlist)  
        input = newsimlist{isim};
        dataFolder    = input.dataFolder;
        dataDirectory = fullfile(mrstOutputDirectory(), directory);
        states = ResultHandler('dataFolder', dataFolder, ...
                               'dataDirectory', dataDirectory);
        
        states = states(:);
        
        ind = cellfun(@(x) not(isempty(x)), states); 
        states = states(ind);
        time = cellfun(@(x) x.time, states); 
        Enew = cellfun(@(x) x.Cathode.E, states); 
        plot((time/hour), Enew, 'linewidth', 3, 'displayname', sprintf('%s, %s %s', input.nparams, input.leaning, directory));
        
    end

end

legend
title('Potential (E)')
xlabel('time (hours)')





%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
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
