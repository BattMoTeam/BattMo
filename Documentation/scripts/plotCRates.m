
close all
figure
hold on

set(0, 'defaulttextfontsize', 15);
set(0, 'defaultaxesfontsize', 15);
set(0, 'defaultlinelinewidth', 3);

CRates = [0.8, 1, 2];
for i = 1 : numel(CRates)
    jsonstruct.Control.CRate = CRates(i);
    output = runBatteryJson(jsonstruct);
    plotResult(output);
end

legend({'CRate = 0.8', ...
        'CRate = 1'  , ...
        'CRate = 2'})

axis([0, 1.5, 2.5, 4.5])


