jsonstruct = parseBattmoJson('exampleECM.json');
inputparams = EquivalentCircuitModelInputParams(jsonstruct);
model = EquivalentCircuitModel(inputparams);

model.solve()
