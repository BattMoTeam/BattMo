function config = configBruggeman(jsonExp)

    config = table();

    % Variable names
    config.name = {
        'ne_am_bruggeman';
        'elyte_bruggeman';
        'sep_bruggeman';
        'pe_am_bruggeman'
                  };
    config.Row = config.name;

    % Variable locations
    numVars = numel(config.name);
    config.belongsTo = repmat({'model'}, numVars, 1);
    config.location = {
        {'NegativeElectrode', 'ActiveMaterial', 'BruggemanCoefficient'};
        {'Electrolyte', 'BruggemanCoefficient'};
        {'Electrolyte', 'Separator', 'BruggemanCoefficient'};
        {'PositiveElectrode', 'ActiveMaterial', 'BruggemanCoefficient'}
                      };

    % Assume Bruggeman coefficients in range
    config.boxLims = repmat([1, 3], numVars, 1);

    % Linear scaling
    config.scaling = repmat({'linear'}, numVars, 1);

    % Store reference values
    pExp = zeros(numVars, 1);
    for k = 1:numVars
        loc = config.location{k};
        pExp(k) = getfield(jsonExp, loc{:});
        config.referenceValue(k) = pExp(k);
    end

    % Initial guess
    p0 = repmat(2, numVars, 1);
    config.initialGuess = p0;

end
