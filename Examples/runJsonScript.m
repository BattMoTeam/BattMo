mrstDebug(20);

mrstModule add ad-core mrst-gui mpfa agmg linearsolvers


% load json struct for the material properties
jsonfilename = fullfile('ParameterData', 'BatteryCellParameters', 'LithiumIonBatteryCell', ...
                        'lithium_ion_battery_nmc_graphite.json');

jsonstruct_material = parseBattmoJson(jsonfilename);

% load json struct for the geometrical properties
% jsonfilename = fullfile('Examples','utils', 'data', 'geometry1d.json');
jsonfilename = './utils/data/4680-geometry.json';
jsonstruct_geometry = parseBattmoJson(jsonfilename);


% merge the two json structs 
jsonstruct = mergeJsonStructs({jsonstruct_geometry,
                               jsonstruct_material}, 'force', true);

% For jellyroll case, we use a iterative linear solver
if strcmp(jsonstruct.Geometry.case, 'jellyRoll')
    
    jsonstruct.NonLinearSolver.maxIterations = 10;
    jsonstruct.NonLinearSolver.verbose = false;

    ne = 'NegativeElectrode';
    pe = 'PositiveElectrode';
    am = 'ActiveMaterial';
    
    diffusionModelType = jsonstruct.(pe).(am).diffusionModelType;
    assert(strcmp(jsonstruct.(pe).(am).diffusionModelType, ...
                  jsonstruct.(ne).(am).diffusionModelType), ...
           'case where two different diffusion models are used in negative and positive electrode has not been implemented.');
    
    switch diffusionModelType
      case 'simple'
        jsonfilename = fullfile(battmoDir, 'Utilities/JsonSchemas/Tests/linearsolver5.json');
      case 'full'
        jsonfilename = fullfile(battmoDir, 'Utilities/JsonSchemas/Tests/linearsolver4.json');
      otherwise
        error('diffusionModelType not covered')
    end

    jsonstruct.NonLinearSolver.LinearSolver.linearSolverSetup = parseBattmoJson(jsonfilename);
    
end

CRate = jsonstruct.Control.CRate;

jsonstruct.TimeStepping.totalTime = 1.4*hour/CRate;
jsonstruct.TimeStepping.N = 100;

% Run battery simulation with function that takes json input
output = runBatteryJson(jsonstruct);

%%

E             = output.E;
energyDensity = output.energyDensity;
energy        = output.energy;

figure
plot(energyDensity, E)
xlabel('Energy Density [Wh/L]');
ylabel('Voltage [V]');

figure
plot(energy, E)
xlabel('Energy [Wh]');
ylabel('Voltage [V]');


