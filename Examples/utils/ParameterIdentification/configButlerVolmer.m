function config = configButlerVolmer(jsonExp)

    ne  = 'NegativeElectrode';
    pe  = 'PositiveElectrode';
    am  = 'ActiveMaterial';
    itf = 'Interface';
    co  = 'Coating';
    k0  = 'reactionRateConstant';

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
    config.location = {{ne, co, am, itf, k0};
                       {pe, co, am, itf, k0}};

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
