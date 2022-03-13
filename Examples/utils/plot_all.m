%% script that plots all the submodels with the potential

%  clean-up states structure
states_new = {};
for i = 1:numel(states)
    if(not(isempty(states{i})))
        states_new{i} = states{i}; 
    end
end
states = states_new; 

%% 
% close all

args = {'field', 'am.phi'};
figure(1), plotResultsSubModel(model, 'pe', states, args{:}), title('pe')
figure(2), plotResultsSubModel(model, 'ne', states, args{:}), title('ne')

argsel = {'field', 'phi'};
figure(3), plotResultsSubModel(model, 'elyte', states, argsel{:}), title('elyte')
figure(4), plotResultsSubModel(model, 'ccpe', states, argsel{:}), title('ccpe')
figure(5), plotResultsSubModel(model, 'ccne', states, argsel{:}), title('ccne')


%{
Copyright 2021-2022 SINTEF Industry, Sustainable Energy Technology
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
