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


