%% Simple plot script used in documentation

close all
figure
hold on

set(0, 'defaulttextfontsize', 15);
set(0, 'defaultaxesfontsize', 15);
set(0, 'defaultlinelinewidth', 3);

DRates = [0.8, 1, 2];
for i = 1 : numel(DRates)
    jsonstruct.Control.DRate = DRates(i);
    output = runBatteryJson(jsonstruct);
    plotResult(output);
end

legend({'DRate = 0.8', ...
        'DRate = 1'  , ...
        'DRate = 2'})

axis([0, 1.5, 2.5, 4.5])


