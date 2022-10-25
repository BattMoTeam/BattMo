function config = configVolumetricSurfaceArea(jsonExp)

    config = table();

    % Variable names
    config.name = {
        'ne_am_volsurfarea';
        'pe_am_volsurfarea'
                  };
    config.Row = config.name;

    % Locations
    numVars = numel(config.name);
    config.belongsTo = repmat({'model'}, numVars, 1);
    config.location = {
        {'NegativeElectrode', 'ActiveMaterial', 'Interface', 'volumetricSurfaceArea'};
        {'PositiveElectrode', 'ActiveMaterial', 'Interface', 'volumetricSurfaceArea'}
                      };

    % Assume surface areas in range
    config.boxLims = repmat([1e5, 1e7], numVars, 1);

    % Assume linear scaling
    config.scaling = repmat({'linear'}, numVars, 1);

    % Store reference values
    pExp = zeros(numVars, 1);
    for k = 1:numVars
        loc = config.location{k};
        pExp(k) = getfield(jsonExp, loc{:});
        config.referenceValue(k) = pExp(k);
    end

    % Initial guess
    p0 = repmat(1e6, numVars, 1);
    config.initialGuess = p0;

end
