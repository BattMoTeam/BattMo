clear all
%close all


set(0, 'defaultlinelinewidth', 2)
set(0, 'defaulttextfontsize', 15);
set(0, 'defaultaxesfontsize', 15);

fd = 'Feed';
pt = 'Permeate';


filename = fullfile(battmoDir(), 'CO2capture/co2membrane.json');
jsonstruct = parseBattmoJson(filename);

inputparams = CO2captureInputParams(jsonstruct);

gen = CO2captureGridGenerator();
gen.nx = 1000;

inputparams = gen.updateInputParams(inputparams);

model = CO2capture(inputparams);



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

channels = {'Feed', 'Permeate'};
grids = {Gf, Gp};
colors = lines(7);

gases = fieldnames(jsonstruct.permeabilities);

% Mol fraction
for igas = 1:numel(gases)
    gas = gases{igas};

    figure, hold on
    for ichannel = 1:2
        channel = channels{ichannel};
        displayname = sprintf('%s %s', channel, gas);

        G = grids{ichannel};
        hp = plot(G.cells.centroids, state.(channel).molFractions{igas}, 'displayname', displayname, 'color', colors(ichannel,:));

        % initstate
        plot(G.cells.centroids, initstate.(channel).molFractions{igas}, 'displayname', sprintf('%s init', displayname), 'color', colors(ichannel,:), 'linestyle', '--');

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
%     for ichannel = 1:2
%         channel = channels{ichannel};
%         displayname = sprintf('%s %s', channel, gas);

%         mf = cellfun(@(x) x.(channel).molFractions{igas}, states, 'uniformoutput', false);
%         mf = cellfun(@(x) sum(x), mf);

%         hp = plot(times, mf, 'displayname', displayname, 'color', colors(ichannel,:));

%         grid on
%     end
%     legend('location', 'sw')
%     xlabel 'time'
%     ylabel 'mol fraction'

% end



% Pressure
figure, hold on
for ichannel = 1:2
    channel = channels{ichannel};
    displayname = sprintf('%s pressure', channel);

    G = grids{ichannel};
    hp = plot(G.cells.centroids, state.(channel).pressure, 'displayname', displayname, 'color', colors(ichannel,:));

    % initstate
    plot(G.cells.centroids, initstate.(channel).pressure, 'displayname', sprintf('%s init', displayname), 'color', colors(ichannel,:), 'linestyle', '--');

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
