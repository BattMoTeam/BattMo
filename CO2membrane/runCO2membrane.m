clear all
%close all


set(0, 'defaultlinelinewidth', 2)
set(0, 'defaulttextfontsize', 15);
set(0, 'defaultaxesfontsize', 15);

fd = 'Feed';
pt = 'Permeate';


filename = fullfile(battmoDir(), 'CO2membrane/co2membrane.json');
jsonstruct = parseBattmoJson(filename);

inputparams = CO2membraneInputParams(jsonstruct);

gen = CO2membraneGridGenerator();
gen.nx = 1000;

inputparams = gen.updateInputParams(inputparams);

model = CO2membrane(inputparams);



shortNames =  {'Feed' ,'fd'  ;
               'Permeate' ,'pt'  ;
               'Boundary' ,'bd'  ;
               'massConses' ,'mcs'  ;
               'bcMolFractionDefinitions', 'bcmfdef';
               'bcFluxDefinition', 'bcfluxdef'};


model = model.equipModelForComputation('shortNames', shortNames);

initstate = model.setupInitialState();

Nt = 100;
dt = 0.1;
step = struct('val', repmat(dt, Nt, 1), ...
              'control', ones(Nt, 1));

control = struct('src', []);

schedule = struct('step', step, ...
                  'control', control);
model.verbose = false;

[~, states] = simulateScheduleAD(initstate, model, schedule);


%% plot
%close all


stateno = numel(states);
state = states{stateno};
state = model.addVariables(state);

Gf = model.Feed.G.mrstFormat();
Gp = model.Permeate.G.mrstFormat();

figure, plotToolbar(Gf, states);

sides = {'Feed', 'Permeate'};
grids = {Gf, Gp};
colors = lines(7);

gases = fieldnames(jsonstruct.permeabilities);

% Mol fraction
for igas = 1:numel(gases)
    gas = gases{igas};

    figure, hold on
    for iside = 1:2
        side = sides{iside};
        displayname = sprintf('%s %s', side, gas);

        G = grids{iside};
        hp = plot(G.cells.centroids, state.(side).molFractions{igas}, 'displayname', displayname, 'color', colors(iside,:));

        % initstate
        plot(G.cells.centroids, initstate.(side).molFractions{igas}, 'displayname', sprintf('%s init', displayname), 'color', colors(iside,:), 'linestyle', '--');

        grid on
    end
    legend('location', 'sw')
    xlabel 'x'
    ylabel 'mol fraction'

end

% % Mol fractions over time
% times = cellfun(@(x) x.time, states);
% for igas = 1:numel(gases)
%     gas = gases{igas};

%     figure, hold on
%     for iside = 1:2
%         side = sides{iside};
%         displayname = sprintf('%s %s', side, gas);

%         mf = cellfun(@(x) x.(side).molFractions{igas}, states, 'uniformoutput', false);
%         mf = cellfun(@(x) sum(x), mf);

%         hp = plot(times, mf, 'displayname', displayname, 'color', colors(iside,:));

%         grid on
%     end
%     legend('location', 'sw')
%     xlabel 'time'
%     ylabel 'mol fraction'

% end



% Pressure
figure, hold on
for iside = 1:2
    side = sides{iside};
    displayname = sprintf('%s pressure', side);

    G = grids{iside};
    hp = plot(G.cells.centroids, state.(side).pressure, 'displayname', displayname, 'color', colors(iside,:));

    % initstate
    plot(G.cells.centroids, initstate.(side).pressure, 'displayname', sprintf('%s init', displayname), 'color', colors(iside,:), 'linestyle', '--');

    grid on
end
legend('location', 'sw')
xlabel 'x'
ylabel 'pressure'


return


cgp = model.cgp
cgt = model.cgt

close all
figure
cgp.plot;

cgt.printRootVariables;
cgt.printTailVariables;
