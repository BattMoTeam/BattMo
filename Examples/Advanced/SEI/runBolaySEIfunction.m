function output = runBolaySEIfunction(input, varargin)
    
    input_default = struct('SOC'  , 1, ...
                           'DiffE', 3.5e-15);

    if ~isempty(input)
        fds = fieldnames(input);
        vals = cellfun(@(fd) input.(fd), fds, 'un', false);
        input = horzcat(fds, vals);
        input = reshape(input', [], 1);
        input = merge_options(input_default, input{:});
    else
        input = input_default;
    end

    ne    = 'NegativeElectrode';
    pe    = 'PositiveElectrode';
    am    = 'ActiveMaterial';
    co    = 'Coating';
    sd    = 'SolidDiffusion';
    itf   = 'Interface';
    sei   = 'SolidElectrodeInterface';
    sr    = 'SideReaction';
    elyte = 'Electrolyte';
    ctrl  = 'Control';

    %% Setup the properties of the Li-ion battery materials and of the cell design
    jsonfilename = fullfile('ParameterData', 'BatteryCellParameters', 'LithiumIonBatteryCell', ...
                            'lithium_ion_battery_nmc_graphite_bolay.json');
    jsonstruct = parseBattmoJson(jsonfilename);

    jsonstruct.use_thermal = false;

    jsonfilename = fullfile('ParameterData', 'ParameterSets', 'Bolay2022', 'bolay_sei_interface.json');
    jsonstruct_bolay = parseBattmoJson(jsonfilename);

    jsonstruct.(ne).(co).(am) = mergeJsonStructs({jsonstruct.(ne).(co).(am), ...
                                                  jsonstruct_bolay});

    jsonstruct.(ne).(co).(am).SEImodel = 'Bolay';

    jsontruct_control = struct( 'controlPolicy'     , 'CCDischarge', ...
                                'initialControl'    , 'discharging', ...
                                'DRate'             , 1            , ...
                                'lowerCutoffVoltage', 3            , ...
                                'rampupTime'        , 100          , ...
                                'upperCutoffVoltage', 4);

    jsonstruct.(ctrl) = jsontruct_control;


    jsonstruct.(ne).(co).(am).(itf).SEIelectronicDiffusionCoefficient = input.DiffE;

    jsonstruct.SOC = input.SOC;

    inputparams = BatteryInputParams(jsonstruct);

    gen = BatteryGeneratorP2D();
    inputparams = gen.updateBatteryInputParams(inputparams);

    model = GenericBattery(inputparams);

    %% Setup the schedule
    %
    model.Control.Imax = 0;

    % schedule = model.(ctrl).setupSchedule([]);

    jsonstruct.TimeStepping.numberOfTimeSteps = 200;
    jsonstruct.TimeStepping.totalTime         = 1*year;

    schedule = model.(ctrl).setupSchedule(jsonstruct);

    %% Setup the initial state of the model
    % The initial state of the model is setup using the model.setupInitialState() method.

    initstate = model.setupInitialState();

    %% Setup non-linear solver

    nls = NonLinearSolver();

    nls.errorOnFailure   = false;
    nls.maxTimestepCuts  = 10;
    nls.timeStepSelector = StateChangeTimeStepSelector('TargetProps', {{'Control','E'}}, 'targetChangeAbs', 0.03);

    model.nonlinearTolerance = 1e-5;

    %% Run simulation

    model.verbose = true;
    [~, states, report] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls);


    %% Setup for plotting

    ind = cellfun(@(x) not(isempty(x)), states); 
    states = states(ind);
    time = cellfun(@(x) x.time, states); 
    E    = cellfun(@(x) x.Control.E, states); 
    I    = cellfun(@(x) x.Control.I, states);

    for istate = 1 : numel(states)
        states{istate} = model.addVariables(states{istate});
    end

    quantities = [];

    vols = model.(ne).(co).G.getVolumes();

    for timeindex = 1 : numel(states)

        state = states{timeindex};
        state = model.evalVarName(state, {ne, co, am, itf, 'SEIconcentration'});
        cSEI  = state.(ne).(co).(am).(itf).SEIconcentration;
        
        Liqqt = sum(cSEI.*vols);
        quantities(end + 1) = Liqqt;
        
    end

    PE_Li_quantities          = [];
    NE_Li_quantities          = [];
    Electrolyte_Li_quantities = [];
    Electrodes_Li_quantities  = [];
    Total_Li_quantities       = [];

    for timeindex = 1 : numel(states)


        amvf     = model.(pe).(co).volumeFractions(1);
        vf       = model.(pe).(co).volumeFraction;
        vols     = model.(pe).G.getVolumes;
        cAverage = states{timeindex}.(pe).(co).(am).(sd).cAverage;

        PE_qtt = sum(amvf.*vf.*vols.*cAverage);
        
        amvf     = model.(ne).(co).volumeFractions(1);
        vf       = model.(ne).(co).volumeFraction;
        vols     = model.(ne).G.getVolumes;
        cAverage = states{timeindex}.(ne).(co).(am).(sd).cAverage;

        NE_qtt = sum(amvf.*vf.*vols.*cAverage);

        Elyte_qtt = sum(model.Electrolyte.volumeFraction.*model.Electrolyte.G.getVolumes.*states{timeindex}.Electrolyte.c);

        Elode_qtt = PE_qtt + NE_qtt;
        Tot_Liqqt = PE_qtt + NE_qtt + Elyte_qtt;
	
        PE_Li_quantities(end + 1)          = PE_qtt;
        NE_Li_quantities(end + 1)          = NE_qtt;
        Electrolyte_Li_quantities(end + 1) = Elyte_qtt;
        Electrodes_Li_quantities(end + 1)  = Elode_qtt;
        Total_Li_quantities(end + 1)       = Tot_Liqqt;

    end

    output = struct('model', model, ...
                    'PE_Li_quantities'         , PE_Li_quantities, ...
                    'NE_Li_quantities'         , NE_Li_quantities, ...
                    'Electrolyte_Li_quantities', Electrolyte_Li_quantities, ...
                    'Electrodes_Li_quantities' , Electrodes_Li_quantities, ...
                    'Total_Li_quantities'      , Total_Li_quantities, ...
                    'quantities'               , quantities, ...
                    'time'                     , time, ...
                    'E'                        , E, ...
                    'I'                        , I);
end
