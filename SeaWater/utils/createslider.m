
if ~exist('dosetup') || dosetup
    dosetup = true;
end

if dosetup
    close all
    states = problem.OutputHandlers.states(:);
    voltagefig = figure(1);
    casenames = {'volumefractions', 'concentrations', 'precipitationrate'};
    mainfigs = {};
    for ind = 1 : numel(casenames)
        mainfigs{end + 1} = figure(ind + 1);
    end
    
    ph = figure('Name'       ,'Time slider', ...
                'NumberTitle','off'        , ...
                'ToolBar'    ,'none'       , ...
                'MenuBar'    ,'none');

    loc = get(ph, 'position');
    loc(4) = 0.5*loc(3);
    set(ph, 'position', loc);
    dosetup = false;
else
    figure(ph)
    clf
end

grp = uibuttongroup('Parent', ph, 'Units', 'Normalized', 'Position', [0 0 1 1]);

t = [1 : numel(states)];

clear input
input.states   = states;
input.model    = model;
input.mainfigs  = mainfigs;
input.casenames = casenames;
input.voltagefig = voltagefig;

timeSlider =  uicontrol('Parent'  , grp                  , ...
                        'Style'   , 'slider'             , ...
                        'min'     , min(t)               , ...
                        'max'     , max(t)               , ...
                        'value'   , 1                    , ...
                        'Units'   , 'normalized'         , ...
                        'Position', [0.1, 0.1, 0.8, 0.5], ...
                        'SliderStep', [1/numel(states), 10/numel(states)], ...
                        'CallBack', {@plotGivenTime, input});

overtext = uicontrol('Parent', grp         , ...
                     'Style' , 'Text'      , ...
                     'Units' , 'Normalized', ...
                     'fontsize', 18, ...
                     'Position', [0.3, 0.8, 0.4, 0.1], ...
                     'String'  , 'time');

set(0, 'defaultlinelinewidth', 3);

% addlistener(timeSlider, 'Value', 'PostSet', @plotGivenTime);




%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
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
