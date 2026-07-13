function [time_sim, I_p2d, V_p2d, ocv_vec] = loadP2dVoltage(time_vec, current_exp)
    
    jsonstruct_material = parseBattmoJson(fullfile('ParameterData','ParameterSets','Chen2020','chen2020_lithium_ion_battery.json'));
    jsonstruct_geometry = parseBattmoJson(fullfile('Examples', 'JsonDataFiles', 'geometryChen.json'));
    jsonstruct = mergeJsonStructs({jsonstruct_material, jsonstruct_geometry});
    
    [model, inputparams] = setupModelFromJson(jsonstruct);
    
    inputparams.Control.controlPolicy = 'CCDischarge'; 
    model.Control = model.setupControl(inputparams.Control);
    
    state_init = setupInitialState(model);
    
    time_vec    = double(time_vec(:));
    current_exp = double(current_exp(:));
    
    dt_steps = diff(time_vec); 
    N_sim_steps = length(dt_steps);
    
    schedule_p2d = struct();
    schedule_p2d.step.val = dt_steps;
    schedule_p2d.step.control = ones(N_sim_steps, 1);
    
    schedule_p2d.control.src = @(t, varargin) interp1(time_vec, current_exp, t, 'linear', 0);
    schedule_p2d.control.Control = struct('stopFunction', @(model, state, state0) false);
    

    [~, states_p2d, ~] = simulateScheduleAD(state_init, model, schedule_p2d);
    
    N_points_sim = length(states_p2d);
    V_p2d    = zeros(N_points_sim, 1);
    I_p2d    = zeros(N_points_sim, 1);
    time_sim = zeros(N_points_sim, 1);
    
    pe_guestStoichiometry0   = model.PositiveElectrode.Coating.ActiveMaterial.Interface.guestStoichiometry0;
    pe_guestStoichiometry100 = model.PositiveElectrode.Coating.ActiveMaterial.Interface.guestStoichiometry100;
    ne_guestStoichiometry0   = model.NegativeElectrode.Coating.ActiveMaterial.Interface.guestStoichiometry0;
    ne_guestStoichiometry100 = model.NegativeElectrode.Coating.ActiveMaterial.Interface.guestStoichiometry100;
    
    for k = 1:N_points_sim
        time_sim(k) = states_p2d{k}.time;
        V_p2d(k)    = states_p2d{k}.Control.E;
        I_p2d(k)    = states_p2d{k}.Control.I;
    end
    
    [Q_nominal, ~] = computeCellCapacity(model); 
    
    soc_p2d = 1 - cumtrapz(time_sim, I_p2d) / Q_nominal;
    soc_p2d = max(0, min(1, soc_p2d));
    
    json_pe_path = fullfile('ParameterData','ParameterSets','Chen2020','chen2020_positive_electrode_interface.json');
    json_pe = parseBattmoJson(json_pe_path);
    [fn_ocp_pe, ~] = setupFunction(json_pe.openCircuitPotential);
    
    json_ne_path = fullfile('ParameterData','ParameterSets','Chen2020','chen2020_negative_electrode_interface.json');
    json_ne = parseBattmoJson(json_ne_path);
    [fn_ocp_ne, ~] = setupFunction(json_ne.openCircuitPotential);
     
    ocv_vec = zeros(N_points_sim, 1);
    for idx = 1:N_points_sim
        soc_now = soc_p2d(idx);
        x_pe = soc_now * (pe_guestStoichiometry100 - pe_guestStoichiometry0) + pe_guestStoichiometry0;              
        x_ne = soc_now * (ne_guestStoichiometry100 - ne_guestStoichiometry0) + ne_guestStoichiometry0;
        ocv_vec(idx) = fn_ocp_pe(x_pe) - fn_ocp_ne(x_ne);
    end
end