function config = configButlerVolmer(jsonExp)

    ne      = 'NegativeElectrode';
    pe      = 'PositiveElectrode';
    elyte   = 'Electrolyte';
    thermal = 'ThermalModel';
    am      = 'ActiveMaterial';
    itf     = 'Interface';
    sd      = 'SolidDiffusion';
    ctrl    = 'Control';
    cc      = 'CurrentCollector';

    config = table();

    % Setup some names for conveneince
    config.name = {'ne_am_k0'; 'pe_am_k0'};
    config.Row = config.name;
    numVars = numel(config.name);

    % The variables are logarithmic
    config.scaling = repmat({'log'}, numVars, 1);

    % Assume the variables belong to some range
    config.boxLims = repmat([1e-12, 1e-9], numVars, 1);

    % The variables belong to the model at a specific location
    config.belongsTo = repmat({'model'}, numVars, 1);
    config.location = {{ne, am, itf, 'k0'}; {pe, am, itf, 'k0'}};

    % Store reference values
    pExp = zeros(numVars, 1);
    for k = 1:numVars
        loc = config.location{k};
        pExp(k) = getfield(jsonExp, loc{:});
        config.referenceValue(k) = pExp(k);
    end

    % Initial guess
    p_0 = [1e-11; 1e-11];
    config.initialGuess = p_0;

end
