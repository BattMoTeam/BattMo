function [Z_real, Z_imag, omegas] = load_chen_data()

        
    % We define some shorthand names for simplicity.
    ne      = 'NegativeElectrode';
    pe      = 'PositiveElectrode';
    elyte   = 'Electrolyte';
    thermal = 'ThermalModel';
    co      = 'Coating';
    am      = 'ActiveMaterial';
    itf     = 'Interface';
    sd      = 'SolidDiffusion';
    ctrl    = 'Control';
    cc      = 'CurrentCollector';
    
    mrstModule add ad-core mrst-gui mpfa agmg linearsolvers
    
    jsonstruct_material = parseBattmoJson(fullfile('ParameterData','ParameterSets','Chen2020','chen2020_lithium_ion_battery.json'));
    jsonstruct_geometry = parseBattmoJson(fullfile('Examples', 'JsonDataFiles', 'geometryChen.json'));
    
    jsonstruct = mergeJsonStructs({jsonstruct_material, ...
                                   jsonstruct_geometry});


    
    includeDoubleLayer = true;
    
    if includeDoubleLayer
    
        jsonstruct.(ne).(co).(am).(itf).useDoubleLayerCapacity = true;
        jsonstruct.(ne).(co).(am).(itf).doubleLayerCapacitance = 0.2;
    
    end
    
    [model, inputparams, ~] = setupModelFromJson(jsonstruct);
    
    c_ne = 29.866*mol/litre; % initial concentration at negative electrode
    c_pe = 17.038*mol/litre; % initial concentration at positive electrode
    
    initstate = initStateChen2020(model, c_ne, c_pe);

    options = [];
    options.stateInitialization.initializationSetup = 'given state';
    options.stateInitialization.computeSteadyState = false;

    extrastructs = [];
    extrastructs.initstate = initstate;
    
    impsolv = ImpedanceSolver(inputparams, options, extrastructs);
    
    omegas = logspace(-4, 2, 50);
    Z = impsolv.computeImpedance(omegas);

    Z_real = real(Z);
    Z_imag = imag(Z);
end
