function efficiencies = computeEfficiency(model, states)
    
    threshold = 1e-5;
    
    ctrl = 'Control';

    E    = cellfun(@(state) state.(ctrl).E, states);
    I    = cellfun(@(state) state.(ctrl).I, states);
    time = cellfun(@(state) state.time, states);

    % ctrlval = 1 if discharging and ctrlval = 2 if charging
    ctrlval = cellfun(@(state) convertCtrlType(state.(ctrl).ctrlType), states); 
    
    assert(isa(model.(ctrl), 'CcCvControlModel'), 'Energy efficiency is only computed for CcCv schedule')

    Umax = model.(ctrl).upperCutoffVoltage;
    Umin = model.(ctrl).lowerCutoffVoltage;

    switch model.(ctrl).initialControl
        
      case 'discharging'

        if (E(1) > Umax)
            
            % We capture the first time when the voltage becomes lower than Umax
            ind = find(E <= Umax, 1);

            % We interpolate the value to match Umax
            c = (Umax - E(ind))/(E(ind - 1) - E(ind));

            X    = assignToX(E, I, time, ctrlval);
            newX = (1 - c)*X(ind, :) + c*X(ind - 1, :);
            X    = [newX; X(ind, :)];

        elseif E(1) < Umax
            
            % Look at the next control
            ind = find(ctrlval == 2, 1);

            % We interpolate to Umin (may be extrapolating)
            c = (Umin - E(ind))/(E(ind - 1) - E(ind));
                
            X    = assignToX(E, I, time, ctrlval);
            newX = (1 - c)*X(ind, :) + c*X(ind - 1, :);
            X    = [newX; X(ind, :)];

        end

      case 'charging'
        
        if (E(1) < Umin)
            % We capture the first time when the voltage becomes larger the Umin
            ind = find(E >= Umin, 1);

            % We interpolate the value to match Umin
            c = (Umin - E(ind))/(E(ind - 1) - E(ind));

            X    = assignToX(E, I, time, ctrlval);
            newX = (1 - c)*X(ind, :) + c*X(ind - 1, :);
            X    = [newX; X(ind, :)];
            
        elseif E(1) > Umin
            
            % Look at the next control
            ind = find(ctrlval == 1, 1);

            % We interpolate to Umin (may be extrapolating)
            c = (Umax - E(ind))/(E(ind - 1) - E(ind));
                
            X    = assignToX(E, I, time, ctrlval);
            newX = (1 - c)*X(ind, :) + c*X(ind - 1, :);
            X    = [newX; X(ind, :)];

        end
        
      otherwise
        
        error('not recognized');
        
    end

    failure = false;
    efficiencies = [];
    
    while ~failure

        output = computeEfficiencyCycle(time, E, I, ctrlval);

        failure = output.failure;

        if failure
            fprintf('%s', output.failureMessage);
            return
        end

        efficiencies(end + 1) = output.efficiencies;

        time    = output.time;
        E       = output.E;
        I       = output.I;
        ctrlval = output.ctrlval;
        
    end
    
end

function output = computeEfficiencyCycle(time, E, I, ctrlval)

    % define threshold to check equality
    threshold = 1e-5;
    
    % We take the first two times where we switch control (from charging to discharging or opposite)
    switches = find(diff(ctrlval) ~= 0, 2, 'first');

    % We check that we obtain a complete cycle

    if isempty(switches) || numel(switches) < 2

        output.isCycleComplete = false;
        output.failure         = true;
        output.failureMessage  = 'Did not find a complete cycle'
        return
        
    end

    if abs(E(1) - E(switches(2) + 1)) > threshold
        
        output.failure = true;
        output.failureMessage = 'The cycle appears to be malformed';
        return
        
    end

    dt   = diff(time);
    Emid = (E(2 : end) + E(1 : end - 1))/2;
    Imid = (I(2 : end) + I(1 : end - 1))/2;

    s1 = 1;
    s2 = switches(1);
    energy1 = sum(Emid(s1 : s2).*Imid(s1 : s2).*dt(s1 : s2));
    
    s1 = switches(1) + 1;
    s2 = switches(2);
    energy2 = sum(Emid(s1 : s2).*Imid(s1 : s2).*dt(s1 : s2));


    if energy1 > 0
        % the first phase is a discharging phase
        efficiency = energy1/abs(energy2);
    else
        efficiency = energy2/abs(energy1);        
    end


    output.efficiency = efficiency;
    output.failure    = false;

    [output.E, output.I, output.time, output.ctrlval] = resetToIndex(switches(2) + 1, E, I, time, ctrlval);
    
end

function [E, I, time, ctrlval] =  resetToIndex(ind, E, I, time, ctrlval, varargin)

    opt = struct('doextrapolation', false)
    opt = merge_options(opt, varargin{:});
    
    E       = E(ind + 1 : end);
    I       = I(ind + 1 : end);
    time    = time(ind + 1 : end);
    ctrlval = ctrlval(ind + 1 : end);
    
end

function ctrlval = convertCtrlType(ctrltype)
    
    switch ctrltype
      case 'CC_discharge1'
        ctrlval = 1;
      case 'CC_discharge2'
        ctrlval = 1;
      case 'CC_charge1'
        ctrlval = 2;
      case 'CV_charge2'
        ctrlval = 2;
      otherwise
        error('ctrlType not recognized');
    end
    
end


function X = assignToX(E, I, time, ctrlval)

    X = [E, I, time, ctrval];
     
end

function [E, I, time, ctrval] = assignFromX(X)
    
    E      = X(:, 1);
    I      = X(:, 2);
    time   = X(:, 3);
    ctrval = X(:, 4);
    
end
