function experience = plotexperience()

    
    json_jp3 = 'C:\Users\Alexandre Fichter\Documents\stage_3A\contenu stage\data_August\jp3-params\jp3-opt-1d-full.json';
    jsonstruct = parseBattmoJson(json_jp3);

    % 1. Construction du modèle
    [model, inputparams] = setupModelFromJson(jsonstruct);
    
    model.NegativeElectrode.Coating.ActiveMaterial.Interface.computeOCPFunc.argumentList = {'c'};
    model.NegativeElectrode.Coating.ActiveMaterial.Interface.computeOCPFunc.numberOfArguments = 1;

    model.PositiveElectrode.Coating.ActiveMaterial.Interface.computeOCPFunc.argumentList = {'c'};
    model.PositiveElectrode.Coating.ActiveMaterial.Interface.computeOCPFunc.numberOfArguments = 1;

    model.Electrolyte.computeDiffusionCoefficientFunc.argumentList = {'c', 'T'};
    model.Electrolyte.computeConductivityFunc.argumentList = {'c', 'T'};
    
    % 3. Initialisation de l'état (maintenant le modèle sait quoi envoyer !)
    state0 = setupInitialState(model);

    inputparams.Control.controlPolicy = 'CCDischarge';
    inputparams.Control.DRate = 2;
    inputparams.Control.tmax = 30;
    inputparams.Control.rampupTime = 0.1;
    inputparams.Control.lowerCutoffVoltage = 2.5;

    schedule1 = setupSchedule(model, inputparams);
    [~, states1] = simulateBattery(model, state0, schedule1);
    state_after_pulse = states1{end};

    inputparams.Control.DRate = 0; 
    inputparams.Control.tmax = 600;
    inputparams.Control.rampupTime = 0.1;


    schedule2 = setupSchedule(model, inputparams);
    disp('Simulation du Repos en cours...');
    [~, states2] = simulateBattery(model, state_after_pulse, schedule2);

    all_states = [states1(:); states2(:)];
    N = length(all_states);
    
    time = zeros(N, 1);
    voltage = zeros(N, 1);
    current = zeros(N, 1);


    for i = 1:N
        time(i) = all_states{i}.time;
        try
            % Dans BattMo, les données du circuit extérieur sont dans "Control"
            voltage(i) = all_states{i}.Control.E;
            current(i) = all_states{i}.Control.I;
        catch
            % Alternative si les variables sont stockées à la racine
            voltage(i) = all_states{i}.voltage;
            current(i) = all_states{i}.current;
        end
    end

    figure;
    subplot(2,1,1);
    plot(time, voltage);
    xlabel('Time (s)');
    ylabel('Voltage (V)');
    title('Voltage vs Time');
    
    subplot(2,1,2);
    plot(time, current);
    xlabel('Time (s)');
    ylabel('Current (A)');
    title('Current vs Time');
    
    experience = [current, voltage, time];
end
