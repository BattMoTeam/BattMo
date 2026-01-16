%% parameter inspection

jsonstruct_mat  = parseBattmoJson('Electrolyser/Parameters/alkalineElectrolyser.json');
jsonstruct_geom = parseBattmoJson('Electrolyser/Parameters/electrolysergeometry1d.json');

% we get a dictionary (struct) with all the parameters
jsonstruct = mergeJsonStructs({jsonstruct_mat, jsonstruct_geom});

% print to screen all the parameter in the dictionary jsonstruct
viewJsonStruct(jsonstruct);

fjv = flattenJsonStruct(jsonstruct);

% use filter option to read subset of parameter data
fjv.print('filter', {'parameter name', 'IonomerMembrane'});

%% setup computation graph interactive tool (cgit)

cgit = model.cgit;

% print all the variables registered in the (computation graph of the) model.
cgit.prinVarNames;

% you can print selections
cgit.printVarNames('Iono');

% write help text in the command window
cgit.help_repl

% given a variable, you can look up the function where it is updated using, for example,
cgit.openPropFunction('Ion conduc');

% plot the graph (in this case no label for the nodes because the graph is too big)
cgit.plot()

% for interactive use
%
figure
cgit.select('OxygenEvolutionElectrode.PorousTransportLayer.eSource')
cgit.parents
cgit.select('PorousT')
cgit.and
cgit.plot

%% processing the results
% they are stored in states

% In runElectrolyser, we have used method addVariables to populate the state for each time-step with all the variables

numel(states) % = 105, one cell per time step, we get here 105 time steps

% we can recover the conductivity for the first time step in the HydrogenEvolutionElectrode.PorousTransportLayer

kappa = state.HydrogenEvolutionElectrode.PorousTransportLayer.conductivity
figure
plot(kappa)

% We get on value per discretization cell

G = model.HydrogenEvolutionElectrode.PorousTransportLayer.grid;
% we can the plot the value of the conductivity at each cell centrid
figure
plot(G.cells.centroids(:,1), state.HydrogenEvolutionElectrode.PorousTransportLayer.conductivity)

%% concentration plot

state = states{end};

% we extract structure that contains index of the different liquid species
liquidInd = model.HydrogenEvolutionElectrode.PorousTransportLayer.liquidInd

figure
hold on

comps = {'H2O', 'OH'};
for icomp = 1 : numel(comps)

    % we plot the concentration of the species comps{icomp}
    c = state.HydrogenEvolutionElectrode.PorousTransportLayer.concentrations{icomp};
    plot(G.cells.centroids(:,1), c, 'displayname', comps{icomp});

end

legend

%% 

% plot the OH concentration at all time steps, on the same plot.
figure
hold on

for istate = 1 : numel(states)
    state = states{istate};
    t = state.time/3600;
    c = state.HydrogenEvolutionElectrode.PorousTransportLayer.concentrations{liquidInd.OH};
    plot(G.cells.centroids(:,1), c, 'displayname', sprintf('%g', t));
end
