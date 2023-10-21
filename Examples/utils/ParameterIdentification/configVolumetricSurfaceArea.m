function config = configVolumetricSurfaceArea(jsonExp)

    ne  = 'NegativeElectrode';
    pe  = 'PositiveElectrode';
    am  = 'ActiveMaterial';
    itf = 'Interface';
    co  = 'Coating';
    vsa = 'volumetricSurfaceArea';

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
    config.location = {{ne, co, am, itf, vsa};
                       {pe, co, am, itf, vsa}};

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



%{
Copyright 2021-2023 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
