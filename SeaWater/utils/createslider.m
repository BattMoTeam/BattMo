
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

