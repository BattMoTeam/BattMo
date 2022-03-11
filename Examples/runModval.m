%%
dorun = false
mrstModule add ad-core mrst-gui mpfa agmg
options = struct();
options.CRate = 1;
options.simcase = 'charge';
options.nas = 20;
options.nL = 10;
options.tabcase = 'aligned tabs'


allopts = {}
tabcases = {'no tab','aligned tabs'}
simcases = {'charge','discharge'}
CRates = [1];
for i=1:numel(simcases)
    for j=1:numel(tabcases)
        for k = 1:numel(CRates)
        opts=options;
        opts.simcase = simcases{i};
        opts.tabcase = tabcases{j};
        opts.CRate = CRates(k);
        allopts{end+1} = opts;
        end
    end
end
figure(1),clf,hold on
if(dorun)
    parfor i = 1:numel(allopts)
        problem = runJellyRollFuncNew(allopts{i});
        %clearPackedSimulatorOutput(problem);
        %simulatePackedProblem(problem)
    end
else
    for i = 1:numel(allopts)
        %problem = runJellyRollFuncNew(allopts{i});
        %[globvars, states, report] = getPackedSimulatorOutput(problem);
        options = allopts{i};
        CRate = options.CRate;
        simcase = options.simcase;
        nas = options.nas;
        nL = options.nL;
        tabcase = options.tabcase;
        filename = ['jellroll_C_',num2str(CRate),'_',simcase,'_nas_',num2str(nas),'_nL_',num2str(nL),'_',tabcase];
        states = ResultHandler('dataPrefix', ['state'], ...
                                          'writeToDisk', true,...
                                          'dataDirectory', '/data/hnil/BITBUCKET/batmo/MRST/mrst-core/output/BatMo/',...
                                          'dataFolder', filename, ...
                                          'cleardir', false);
    
        
    
        states = states(:)
        ind = cellfun(@(x) not(isempty(x)), states); 
        states = states(ind);
        pe = 'PositiveElectrode';
        cc = 'CurrentCollector';
        Enew = cellfun(@(x) x.(pe).(cc).E, states); 
        Inew = cellfun(@(x) x.(pe).(cc).I, states);
        Tmax = cellfun(@(x) max(x.ThermalModel.T), states);
        %[SOCN,SOCP] =  cellfun(@(x) model.calculateSOC(x), states);
        time = cellfun(@(x) x.time, states);
        plot(time,Enew)
    end
end
%%
allproblem={};
allstates={}
for i=1:numel(allopts)
    allproblem{i} = runJellyRollFuncNew(allopts{i});
    [globvars, allstates{i}, report] = getPackedSimulatorOutput(allproblem{i});
end
%%
figure(1),clf
i=2
plotToolbar(allproblem{i}.SimulatorSetup.model.G,allstates{i})
view([0 0 -1])
%%
%%
i =2;
problem2 = runJellyRollFuncNew(allopts{i});
[globvars, states2, report] = getPackedSimulatorOutput(problem2);
%%
figure(2),clf
G2 = problem.SimulatorSetup.model.G;
plotToolbar(G2,states2)
%%
return
%%
figure,hold on
dorun = false;
allopts = {}
tabcases = {'no tab','aligned tabs'}
simcases = {'charge'}
CRates = [0.1];
for i=1:numel(simcases)
    for j=1:numel(tabcases)
        for k = 1:numel(CRates)
        opts=options;
        opts.simcase = simcases{i};
        opts.tabcase = tabcases{j};
        opts.CRate = CRates(k);
        allopts{end+1} = opts;
        end
    end
end
if(dorun)
    parfor i = 1:numel(allopts)
        problem = runJellyRollFuncNew(allopts{i});
        %    clearPackedSimulatorOutput(problem);
        simulatePackedProblem(problem)
    end
else
    for i = 1:numel(allopts)
        %[globvars, states, report] = getPackedSimulatorOutput(problem);
        options = allopts{i};
        CRate = options.CRate;
        simcase = options.simcase;
        nas = options.nas;
        nL = options.nL;
        tabcase = options.tabcase;
        filename = ['jellroll_C_',num2str(CRate),'_',simcase,'_nas_',num2str(nas),'_nL_',num2str(nL),'_',tabcase];
        states = ResultHandler('dataPrefix', ['state'], ...
                                          'writeToDisk', true,...
                                          'dataDirectory', '/data/hnil/BITBUCKET/batmo/MRST/mrst-core/output/BatMo/',...
                                          'dataFolder', filename, ...
                                          'cleardir', false);
    
        
    
        states = states(:)
        ind = cellfun(@(x) not(isempty(x)), states); 
        states = states(ind);
        pe = 'PositiveElectrode';
        cc = 'CurrentCollector';
        Enew = cellfun(@(x) x.(pe).(cc).E, states); 
        Inew = cellfun(@(x) x.(pe).(cc).I, states);
        Tmax = cellfun(@(x) max(x.ThermalModel.T), states);
        %[SOCN,SOCP] =  cellfun(@(x) model.calculateSOC(x), states);
        time = cellfun(@(x) x.time, states);
        plot(time,Enew)
    end
end