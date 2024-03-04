close all
mrstDebug(20);

jsonstruct      = parseBattmoJson('Examples/LNMO/SiGrLnmo.json');
jsonstruct_cccv = parseBattmoJson('Examples/JsonDataFiles/cccv_control.json');

jsonstruct = removeJsonStructFields(jsonstruct, ...
                                    {'Control', 'DRate'}        , ...
                                    {'Control', 'controlPolicy'}, ...
                                    {'Control', 'dEdtLimit'}    , ...
                                    {'Control', 'dIdtLimit'}    , ...
                                    {'Control', 'rampupTime'}   , ...
                                    {'Control', 'useCVswitch'});

jsonstruct = mergeJsonStructs({jsonstruct, jsonstruct_cccv});


ne  = 'NegativeElectrode';
pe  = 'PositiveElectrode';
co  = 'Coating';
am1 = 'ActiveMaterial1';
am2 = 'ActiveMaterial2';
am  = 'ActiveMaterial';
itf = 'Interface';

jsonstruct.(ne).(co).(am1).(itf).guestStoichiometry100 = 1;
jsonstruct.(ne).(co).(am1).(itf).guestStoichiometry0   = 0.01;
jsonstruct.(ne).(co).(am2).(itf).guestStoichiometry100 = 1;
jsonstruct.(ne).(co).(am2).(itf).guestStoichiometry0   = 0.01;

jsonstruct.(pe).(co).(am).(itf).guestStoichiometry0 = 1;

model = setupModelFromJson(jsonstruct);

ecs = EquilibriumConcentrationSolver(model);

[E, extras] = advancedComputeCellEnergy(model);
m = computeCellMass(model);

fprintf('Specific Energy %g\n', (E/m)/hour);

dsfunc = extras.dischargeFunction;
soc = linspace(0, 1, 100);

figure
plot(soc, dsfunc(soc));
xlabel('SOC / -');
ylabel('Voltage / V');

