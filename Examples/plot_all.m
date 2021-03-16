%% script that plots all the submodels with the potential

%  clean-up states structure
states_new = {}
for i = 1:numel(states)
    if(not(isempty(states{i})))
        states_new{i} = states{i}; 
    end
end
states = states_new; 

%% 
% close all

args = {'field', 'am.phi'}
figure(1), plotResultsSubModel(model, 'pe', states, args{:}), title('pe')
figure(2), plotResultsSubModel(model, 'ne', states, args{:}), title('ne')

argsel = {'field', 'phi'}
figure(3), plotResultsSubModel(model, 'elyte', states, argsel{:}), title('elyte')
figure(4), plotResultsSubModel(model, 'ccpe', states, argsel{:}), title('ccpe')
figure(5), plotResultsSubModel(model, 'ccne', states, argsel{:}), title('ccne')