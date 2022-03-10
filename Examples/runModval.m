%%
options = struct();
options.CRate = 1;
options.simcase = 'charge';
options.nas = 5;
options.nL = 2;
options.tabcase = 'aligned tabs'


allopts = {}
tabcases = {'no tab','aligned tabs'}
simcases = {'charge','discharge'}
for i=1:numel(simcases)
    for j=1:numel(tabcases)
        opts=options;
        opts.simcase = simcases{i};
        opts.tabcase = tabcases{j};
        allopts{end+1} = opts;
    end
end

for i = 1:1%numel(allopts)
    problem = runJellyRollFuncNew(allopts{1});
    clearPackedSimulatorOutput(problem);
    simulatePackedProblem(problem)
end