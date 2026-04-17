jsonstruct = parseBattmoJson('exampleECM.json');
inputparams = EquivalentCircuitModelInputParams(jsonstruct);
model = EquivalentCircuitModel(inputparams);

cgit = model.cgit;

clear state;
state.I = 3;

state = model.evalVarName(state, {'dSOCdt'});
state = model.evalVarName(state, {'U'});
