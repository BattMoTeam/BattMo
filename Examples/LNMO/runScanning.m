sc = Scanning();

thickness = linspace(50*micro*meter, 100*micro*meter, 5);
SiContent  = linspace(0.01, 0.1, 5);

[inds, inputs] = cartesianProduct(thickness, SiContent);

clear output

for iinput = 1 : numel(inputs)

    input = inputs{iinput};
    sc.computeSpecificEnergy(input.thickness, input.SiContent);

    E(iinput) = sc.specificEnergyHour;

end

%%

close all

set(0, 'defaultlinelinewidth', 3);
set(0, 'defaultaxesfontsize', 15);

t = cellfun(@(input) input.thickness, inputs);
s = cellfun(@(input) input.SiContent, inputs);

figure
plot3(t, s, E);

X = repmat(thickness', 1, numel(SiContent))/(micro*meter);
Y = repmat(SiContent, numel(thickness), 1);
Z = reshape(E, numel(thickness), numel(SiContent));

figure
surf(X, Y, Z);
xlabel('Negative Electrode Length / μm');
ylabel('Silicon Ratio (mass fraction)');
zlabel('Specific Energy / Wh/kg');
