jsonstruct = parseBattmoJson('exampleECM.json');
inputparams = EquivalentCircuitModelInputParams(jsonstruct);
model = EquivalentCircuitModel(inputparams);

[t, U, I] = model.solve();

%%

tiledlayout(1, 2);
nexttile
plot(t, U)
nexttile
plot(t, I)
