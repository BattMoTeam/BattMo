function [fig] = plotDashboard(model, states, varargin)
%PLOTDASHBOARD Displays a summary of the results for a given simulation.
%   fig = plotDashboard(model, states)
%   
%   fig = plotDashboard(model, states, 'PropertyName1', PropertyValue1)

%
    close all

    %% Parse inputs            
    defaultStep         = length(states);
    defaultTheme        = 'light';
    expectedTheme       = {'dark', 'light', 'blue', 'offwhite'};
    defaultSize         = 'wide';
    expectedSize        = {'A4', 'square', 'wide'};
    defaultOrientation  = 'landscape';
    expectedOrientation = {'landscape', 'portrait'};

    p = inputParser;
    validStep = @(x) isnumeric(x) && isinteger(int8(x)) && (x >= 0) && (x <= length(states));
    addParameter(p, 'step', defaultStep, validStep);
    addParameter(p, 'theme', defaultTheme, ...
                @(x) any(validatestring(x,expectedTheme)));
    addParameter(p, 'size', defaultSize, ...
                @(x) any(validatestring(x,expectedSize)));
    addParameter(p, 'orientation', defaultOrientation, ...
                @(x) any(validatestring(x,expectedOrientation)));
    parse(p, varargin{:});

    % shorthands
    ne      = 'NegativeElectrode';
    pe      = 'PositiveElectrode';
    am      = 'ActiveMaterial';
    cc      = 'CurrentCollector';
    elyte   = 'Electrolyte';
    thermal = 'ThermalModel';
    ctrl    = 'Control';
    sd      = 'SolidDiffusion';
    
    step = p.Results.step;

    fig = figure;

    time = cellfun(@(state) state.time, states); 
    Enew = cellfun(@(state) state.(ctrl).E, states); 
    Inew = cellfun(@(state) state.(ctrl).I, states);


    if step ~= 0
        
        timeBar = [  time(step)/hour, 0; ...
                     time(step)/hour, 1000];
        
        if model.G.griddim == 1
            setFigureStyle('theme', p.Results.theme, 'size', p.Results.size, 'orientation', p.Results.orientation, 'quantity', 'single');
            switch model.(ne).(am).diffusionModelType
              case 'simple'
                subplot(2,4,1), plotCellData(model.(ne).(am).G, states{step}.(ne).(am).c ./ 1000, 'linewidth', 3);
              case 'full'
                subplot(2,4,1), plotCellData(model.(ne).(am).G, states{step}.(ne).(am).(sd).cSurface ./ 1000, 'linewidth', 3);
              otherwise
                error('diffusionModelType not recognized');
            end
            xlabel(gca, 'Position  /  m')
            title(gca, 'Negative Electrode Concentration  /  mol \cdot L^{-1}')
            set(gca, ...
                'FontSize' , style.fontSize       , ...
                'FontName' , style.fontName       , ...
                'color'    , style.backgroundColor, ...
                'XColor'   , style.fontColor      , ...
                'YColor'   , style.fontColor      , ...
                'GridColor', style.fontColor)
            
            subplot(2,4,2), plotCellData(model.(elyte).G, states{step}.(elyte).c ./ 1000, 'linewidth', 3);
            xlabel(gca, 'Position  /  m')
            title(gca, 'Electrolyte Concentration  /  mol \cdot L^{-1}')
            set(gca, ...
                'FontSize' , style.fontSize       , ...
                'FontName' , style.fontName       , ...
                'color'    , style.backgroundColor, ...
                'XColor'   , style.fontColor      , ...
                'YColor'   , style.fontColor      , ...
                'GridColor', style.fontColor)
            
            switch model.(pe).(am).diffusionModelType
              case 'simple'
                subplot(2,4,3), plotCellData(model.(pe).(am).G, states{step}.(pe).(am).c ./ 1000, 'linewidth', 3);
              case 'full'
                subplot(2,4,3), plotCellData(model.(pe).(am).G, states{step}.(pe).(am).(sd).cSurface ./ 1000, 'linewidth', 3);
              otherwise
                error('diffusionModelType not recognized');
            end
            
            xlabel(gca, 'Position  /  m')
            title(gca, 'Positive Electrode Concentration  /  mol \cdot L^{-1}')
            set(gca, ...
                'FontSize' , style.fontSize       , ...
                'FontName' , style.fontName       , ...
                'color'    , style.backgroundColor, ...
                'XColor'   , style.fontColor      , ...
                'YColor'   , style.fontColor      , ...
                'GridColor', style.fontColor)
            
            subplot(2,4,4), plot((time/hour), Inew, '-', 'linewidth', 3)
            hold on
            plot(timeBar(:,1), timeBar(:,2), 'k--', 'linewidth', 1);
            hold off
            title('Cell Current  /  A')
            xlabel('Time  /  h')
            xlim([min(time/hour), max(time/hour)]);
            ylim([min(Inew), max(Inew)]);
            set(gca, ...
                'FontSize' , style.fontSize       , ...
                'FontName' , style.fontName       , ...
                'color'    , style.backgroundColor, ...
                'XColor'   , style.fontColor      , ...
                'YColor'   , style.fontColor      , ...
                'GridColor', style.fontColor)
            
            subplot(2,4,5), plotCellData(model.(ne).(am).G, states{step}.(ne).(am).phi, 'linewidth', 3);
            xlabel(gca, 'Position  /  m')
            title(gca, 'Negative Electrode Potential  /  V')
            set(gca, ...
                'FontSize' , style.fontSize       , ...
                'FontName' , style.fontName       , ...
                'color'    , style.backgroundColor, ...
                'XColor'   , style.fontColor      , ...
                'YColor'   , style.fontColor      , ...
                'GridColor', style.fontColor)
            
            subplot(2,4,6), plotCellData(model.(elyte).G, states{step}.(elyte).phi, 'linewidth', 3);
            xlabel(gca, 'Position  /  m')
            title(gca, 'Electrolyte Potential  /  V')
            set(gca, ...
                'FontSize' , style.fontSize       , ...
                'FontName' , style.fontName       , ...
                'color'    , style.backgroundColor, ...
                'XColor'   , style.fontColor      , ...
                'YColor'   , style.fontColor      , ...
                'GridColor', style.fontColor)
            
            subplot(2,4,7), plotCellData(model.(pe).(am).G, states{step}.(pe).(am).phi, 'linewidth', 3);
            xlabel(gca, 'Position  /  m')
            title(gca, 'Positive Electrode Potential  /  V')
            set(gca, ...
                'FontSize' , style.fontSize       , ...
                'FontName' , style.fontName       , ...
                'color'    , style.backgroundColor, ...
                'XColor'   , style.fontColor      , ...
                'YColor'   , style.fontColor      , ...
                'GridColor', style.fontColor)
            
            subplot(2,4,8), plot((time/hour), Enew, '-', 'linewidth', 3)
            hold on
            plot(timeBar(:,1), timeBar(:,2), 'k--', 'linewidth', 1);
            hold off
            title('Cell Voltage  /  V')
            xlabel('Time  /  h')
            xlim([min(time/hour), max(time/hour)]);
            ylim([min(Enew), max(Enew)]);
            set(gca, ...
                'FontSize' , style.fontSize       , ...
                'FontName' , style.fontName       , ...
                'color'    , style.backgroundColor, ...
                'XColor'   , style.fontColor      , ...
                'YColor'   , style.fontColor      , ...
                'GridColor', style.fontColor)

        else
            style = setFigureStyle('theme', p.Results.theme, 'size', p.Results.size, 'orientation', p.Results.orientation, 'quantity', 'single');
            style.fontSize = 10;

            switch model.(ne).(am).diffusionModelType
              case 'simple'
                subplot(2,4,1), plotCellData(model.(ne).(am).G, states{step}.(ne).(am).c ./ 1000, 'edgealpha', 0.1);
              case 'full'
                subplot(2,4,1), plotCellData(model.(ne).(am).G, states{step}.(ne).(am).(sd).cSurface ./ 1000, 'edgealpha', 0.1);
              otherwise
                error('diffusionModelType not recognized');
            end
            
            xlabel(gca, 'Position  /  m')
            ylabel(gca, 'Position  /  m')
            title(gca, 'Negative Electrode Concentration  /  mol \cdot L^{-1}')
            colormap(crameri('nuuk'));
            colorbar
            axis tight
            set(gca, ...
                'FontSize' , style.fontSize       , ...
                'FontName' , style.fontName       , ...
                'color'    , style.backgroundColor, ...
                'XColor'   , style.fontColor      , ...
                'YColor'   , style.fontColor      , ...
                'GridColor', style.fontColor)
            if model.G.griddim == 3
                view(45,45);
                axis equal
            end
            
            
            
            subplot(2,4,2), plotCellData(model.(elyte).G, states{step}.(elyte).c ./ 1000, 'edgealpha', 0.1);
            xlabel(gca, 'Position  /  m')
            title(gca, 'Electrolyte Concentration  /  mol \cdot L^{-1}')
            colormap(crameri('nuuk'));
            colorbar
            axis tight
            set(gca, ...
                'FontSize' , style.fontSize       , ...
                'FontName' , style.fontName       , ...
                'color'    , style.backgroundColor, ...
                'XColor'   , style.fontColor      , ...
                'YColor'   , style.fontColor      , ...
                'GridColor', style.fontColor)
            if model.G.griddim == 3
                view(45,45);
                axis equal
            end

            switch model.(pe).(am).diffusionModelType
              case 'simple'
                subplot(2,4,3), plotCellData(model.(pe).(am).G, states{step}.(pe).(am).c ./ 1000, 'edgealpha', 0.1);
              case 'full'
                subplot(2,4,3), plotCellData(model.(pe).(am).G, states{step}.(pe).(am).(sd).c ./ 1000, 'edgealpha', 0.1);
              otherwise
                error('diffusionModelType not recognized');
            end
                
            xlabel(gca, 'Position  /  m')
            title(gca, 'Positive Electrode Concentration  /  mol \cdot L^{-1}')
            colormap(crameri('nuuk'));
            colorbar
            axis tight
            set(gca, ...
                'FontSize' , style.fontSize       , ...
                'FontName' , style.fontName       , ...
                'color'    , style.backgroundColor, ...
                'XColor'   , style.fontColor      , ...
                'YColor'   , style.fontColor      , ...
                'GridColor', style.fontColor)
            if model.G.griddim == 3
                view(45,45);
                axis equal
            end
            
            subplot(2,4,4), plot((time/hour), Inew, '-', 'linewidth', 3)
            hold on
            plot(timeBar(:,1), timeBar(:,2), 'k--', 'linewidth', 1);
            hold off
            title('Cell Current  /  A')
            xlabel('Time  /  h')
            xlim([min(time/hour), max(time/hour)]);
            ylim([min(Inew), max(Inew)]);
            set(gca, ...
                'FontSize' , style.fontSize       , ...
                'FontName' , style.fontName       , ...
                'color'    , style.backgroundColor, ...
                'XColor'   , style.fontColor      , ...
                'YColor'   , style.fontColor      , ...
                'GridColor', style.fontColor)
            
            subplot(2,4,5), plotCellData(model.(ne).(am).G, states{step}.(ne).(am).phi, 'edgealpha', 0.1);
            xlabel(gca, 'Position  /  m')
            title(gca, 'Negative Electrode Potential  /  V')
            colormap(gca, crameri('lapaz'))
            colorbar
            axis tight
            set(gca, ...
                'FontSize' , style.fontSize       , ...
                'FontName' , style.fontName       , ...
                'color'    , style.backgroundColor, ...
                'XColor'   , style.fontColor      , ...
                'YColor'   , style.fontColor      , ...
                'GridColor', style.fontColor)
            if model.G.griddim == 3
                view(45,45);
                axis equal
            end
            
            subplot(2,4,6), plotCellData(model.(elyte).G, states{step}.(elyte).phi, 'edgealpha', 0.1);
            xlabel(gca, 'Position  /  m')
            title(gca, 'Electrolyte Potential  /  V')
            colormap(gca, crameri('lapaz'))
            colorbar
            axis tight
            set(gca, ...
                'FontSize' , style.fontSize       , ...
                'FontName' , style.fontName       , ...
                'color'    , style.backgroundColor, ...
                'XColor'   , style.fontColor      , ...
                'YColor'   , style.fontColor      , ...
                'GridColor', style.fontColor)
            if model.G.griddim == 3
                view(45,45);
                axis equal
            end
            
            subplot(2,4,7), plotCellData(model.(pe).(am).G, states{step}.(pe).(am).phi, 'edgealpha', 0.1);
            xlabel(gca, 'Position  /  m')
            title(gca, 'Positive Electrode Potential  /  V')
            colormap(gca, crameri('lapaz'))
            colorbar
            axis tight
            
            set(gca, ...
                'FontSize' , style.fontSize       , ...
                'FontName' , style.fontName       , ...
                'color'    , style.backgroundColor, ...
                'XColor'   , style.fontColor      , ...
                'YColor'   , style.fontColor      , ...
                'GridColor', style.fontColor)
            if model.G.griddim == 3
                view(45,45);
                axis equal
            end
            
            subplot(2,4,8), plot((time/hour), Enew, '-', 'linewidth', 3)
            hold on
            plot(timeBar(:,1), timeBar(:,2), 'k--', 'linewidth', 1);
            hold off
            title('Cell Voltage  /  V')
            xlabel('Time  /  h')
            xlim([min(time/hour), max(time/hour)]);
            ylim([min(Enew), max(Enew)]);
            set(gca, ...
                'FontSize' , style.fontSize       , ...
                'FontName' , style.fontName       , ...
                'color'    , style.backgroundColor, ...
                'XColor'   , style.fontColor      , ...
                'YColor'   , style.fontColor      , ...
                'GridColor', style.fontColor)

        end
        %     setFigureStyle('theme', p.Results.theme, 'size', p.Results.size, 'orientation', p.Results.orientation, 'quantity', 'single');
    else
        
        for i = 1 : length(states)
            
            if i == 1
                
                cmax_elyte = max(max(states{i}.(elyte).c ./ 1000));
                cmin_elyte = min(min(states{i}.(elyte).c ./ 1000));
                
                switch model.(ne).(am).diffusionModelType
                  case 'simple'
                    cmax_ne = max(max(states{i}.(ne).(am).c ./ 1000));
                    cmin_ne = min(min(states{i}.(ne).(am).c ./ 1000));
                  case 'full'
                    cmax_ne = max(max(states{i}.(ne).(am).(sd).cSurface ./ 1000));
                    cmin_ne = min(min(states{i}.(ne).(am).(sd).cSurface ./ 1000));
                  otherwise
                    error('diffusionModelType not recognized');
                end

                switch model.(pe).(am).diffusionModelType
                  case 'simple'
                    cmax_pe = max(max(states{i}.(pe).(am).c ./ 1000));
                    cmin_pe = min(min(states{i}.(pe).(am).c ./ 1000));
                  case 'full'
                    cmax_pe = max(max(states{i}.(pe).(am).(sd).cSurface ./ 1000));
                    cmin_pe = min(min(states{i}.(pe).(am).(sd).cSurface ./ 1000));
                  otherwise
                    error('diffusionModelType not recognized');
                end
                
                
                phimax_elyte = max(max(states{i}.(elyte).phi));
                phimin_elyte = min(min(states{i}.(elyte).phi));
                
                phimax_ne = max(max(states{i}.(ne).(am).phi));
                phimin_ne = min(min(states{i}.(ne).(am).phi));
                
                phimax_pe = max(max(states{i}.(pe).(am).phi));
                phimin_pe = min(min(states{i}.(pe).(am).phi));
                
                xmin = min(model.(elyte).G.nodes.coords(:,1));
                xmax = max(model.(elyte).G.nodes.coords(:,1));
                if model.G.griddim == 2
                    ymin = min(model.(elyte).G.nodes.coords(:,2));
                    ymax = max(model.(elyte).G.nodes.coords(:,2));
                elseif model.G.griddim == 3
                    ymin = min(model.(elyte).G.nodes.coords(:,2));
                    ymax = max(model.(elyte).G.nodes.coords(:,2));
                    zmin = min(model.(elyte).G.nodes.coords(:,3));
                    zmax = max(model.(elyte).G.nodes.coords(:,3));
                end
            else
                
                cmax_elyte = max(cmax_elyte, max(max(states{i}.(elyte).c ./ 1000)));
                cmin_elyte = min(cmin_elyte, min(min(states{i}.(elyte).c ./ 1000)));
                
                switch model.(ne).(am).diffusionModelType
                  case 'simple'
                    cmax_ne = max(cmax_ne, max(max(states{i}.(ne).(am).c ./ 1000)));
                    cmin_ne = min(cmin_ne, min(min(states{i}.(ne).(am).c ./ 1000)));
                  case 'full'
                    cmax_ne = max(cmax_ne, max(max(states{i}.(ne).(am).(sd).cSurface./ 1000)));
                    cmin_ne = min(cmin_ne, min(min(states{i}.(ne).(am).(sd).cSurface./ 1000)));
                  otherwise
                    error('diffusionModelType not recognized');
                end
                switch model.(pe).(am).diffusionModelType
                  case 'simple'
                    cmax_pe = max(cmax_pe, max(max(states{i}.(pe).(am).c ./ 1000)));
                    cmin_pe = min(cmin_pe, min(min(states{i}.(pe).(am).c ./ 1000)));
                  case 'full'
                    cmax_pe = max(cmax_pe, max(max(states{i}.(pe).(am).(sd).cSurface./ 1000)));
                    cmin_pe = min(cmin_pe, min(min(states{i}.(pe).(am).(sd).cSurface./ 1000)));
                  otherwise
                    error('diffusionModelType not recognized');
                end

                cmax_global_solid = max(cmax_ne, cmax_pe);
                cmin_global_solid = min(cmin_ne, cmin_pe);
                
                phimax_elyte = max(phimax_elyte, max(max(states{i}.(elyte).phi)));
                phimin_elyte = min(phimin_elyte, min(min(states{i}.(elyte).phi)));
                
                phimax_ne = max(phimax_ne, max(max(states{i}.(ne).(am).phi)));
                phimin_ne = min(phimin_ne, min(min(states{i}.(ne).(am).phi)));
                
                phimax_pe = max(phimax_pe, max(max(states{i}.(pe).(am).phi)));
                phimin_pe = min(phimin_pe, min(min(states{i}.(pe).(am).phi)));
                
                phimax_global = max([phimax_ne, phimax_pe, phimax_elyte]);
                phimin_global = min([cmin_ne, cmin_pe, phimin_elyte]);
            end
        end
        
        for i = 1:length(states)
            
            timeBar = [  time(i)/hour, 0; ...
                         time(i)/hour, 1000];
            
            if i == 1
                style = setFigureStyle('theme', p.Results.theme, 'size', p.Results.size, 'orientation', p.Results.orientation, 'quantity', 'single'); 
                style.fontSize = 10;
            end
            if model.G.griddim == 1

                switch model.(ne).(am).diffusionModelType
                  case 'simple'
                    subplot(2,4,1), plotCellData(model.(ne).(am).G, states{i}.(ne).(am).c ./ 1000, 'linewidth', 3);
                  case 'full'
                    subplot(2,4,1), plotCellData(model.(ne).(am).G, states{i}.(ne).(am).(sd).cSurface ./ 1000, 'linewidth', 3);
                  otherwise
                    error('diffusionModelType not recognized');
                end
                
                xlabel(gca, 'Position  /  m')
                title(gca, 'Negative Electrode Concentration  /  mol \cdot L^{-1}', 'color', style.fontColor)
                xlim([xmin, xmax])
                ylim([cmin_ne, cmax_ne])
                set(gca, ...
                    'FontSize' , style.fontSize       , ...
                    'FontName' , style.fontName       , ...
                    'color'    , style.backgroundColor, ...
                    'XColor'   , style.fontColor      , ...
                    'YColor'   , style.fontColor      , ...
                    'GridColor', style.fontColor)

                subplot(2,4,2), plotCellData(model.(elyte).G, states{i}.(elyte).c ./ 1000, 'linewidth', 3);

                    
                xlabel(gca, 'Position  /  m')
                title(gca, 'Electrolyte Concentration  /  mol \cdot L^{-1}', 'color', style.fontColor)
                xlim([xmin, xmax])
                ylim([cmin_elyte, cmax_elyte])
                set(gca, ...
                    'FontSize' , style.fontSize       , ...
                    'FontName' , style.fontName       , ...
                    'color'    , style.backgroundColor, ...
                    'XColor'   , style.fontColor      , ...
                    'YColor'   , style.fontColor      , ...
                    'GridColor', style.fontColor)


                switch model.(pe).(am).diffusionModelType
                  case 'simple'
                    subplot(2,4,3), plotCellData(model.(pe).(am).G, states{i}.(pe).(am).c ./ 1000, 'linewidth', 3);
                  case 'full'
                    subplot(2,4,3), plotCellData(model.(pe).(am).G, states{i}.(pe).(am).(sd).cSurface ./ 1000, 'linewidth', 3);
                  otherwise
                    error('diffusionModelType not recognized');
                end
                
                xlabel(gca, 'Position  /  m')
                title(gca, 'Positive Electrode Concentration  /  mol \cdot L^{-1}', 'color', style.fontColor)
                xlim([xmin, xmax])
                ylim([cmin_pe, cmax_pe])
                set(gca, ...
                    'FontSize' , style.fontSize       , ...
                    'FontName' , style.fontName       , ...
                    'color'    , style.backgroundColor, ...
                    'XColor'   , style.fontColor      , ...
                    'YColor'   , style.fontColor      , ...
                    'GridColor', style.fontColor)

                subplot(2,4,4), plot((time/hour), Inew, '-', 'linewidth', 3)
                hold on
                plot(timeBar(:,1), timeBar(:,2),  '--', 'linewidth', 1, 'color', style.fontColor);
                hold off
                title('Cell Current  /  A', 'color', style.fontColor)
                xlabel('Time  /  h')
                xlim([min(time/hour), max(time/hour)]);
                ylim([min(Inew), max(Inew)]);
                set(gca, ...
                    'FontSize' , style.fontSize       , ...
                    'FontName' , style.fontName       , ...
                    'color'    , style.backgroundColor, ...
                    'XColor'   , style.fontColor      , ...
                    'YColor'   , style.fontColor      , ...
                    'GridColor', style.fontColor)

                subplot(2,4,5), plotCellData(model.(ne).(am).G, states{i}.(ne).(am).phi, 'linewidth', 3);
                xlabel(gca, 'Position  /  m')
                title(gca, 'Negative Electrode Potential  /  V', 'color', style.fontColor)
                xlim([xmin, xmax])
                ylim([phimin_ne, phimax_ne])
                set(gca, ...
                    'FontSize' , style.fontSize       , ...
                    'FontName' , style.fontName       , ...
                    'color'    , style.backgroundColor, ...
                    'XColor'   , style.fontColor      , ...
                    'YColor'   , style.fontColor      , ...
                    'GridColor', style.fontColor)

                subplot(2,4,6), plotCellData(model.(elyte).G, states{i}.(elyte).phi, 'linewidth', 3);
                xlabel(gca, 'Position  /  m')
                title(gca, 'Electrolyte Potential  /  V', 'color', style.fontColor)
                xlim([xmin, xmax])
                ylim([phimin_elyte, phimax_elyte])
                set(gca, ...
                    'FontSize' , style.fontSize       , ...
                    'FontName' , style.fontName       , ...
                    'color'    , style.backgroundColor, ...
                    'XColor'   , style.fontColor      , ...
                    'YColor'   , style.fontColor      , ...
                    'GridColor', style.fontColor)

                subplot(2,4,7), plotCellData(model.(pe).(am).G, states{i}.(pe).(am).phi, 'linewidth', 3);
                xlabel(gca, 'Position  /  m')
                title(gca, 'Positive Electrode Potential  /  V', 'color', style.fontColor)
                xlim([xmin, xmax])
                ylim([phimin_pe, phimax_pe])
                set(gca, ...
                    'FontSize' , style.fontSize       , ...
                    'FontName' , style.fontName       , ...
                    'color'    , style.backgroundColor, ...
                    'XColor'   , style.fontColor      , ...
                    'YColor'   , style.fontColor      , ...
                    'GridColor', style.fontColor)

                subplot(2,4,8), plot((time/hour), Enew, '-', 'linewidth', 3)
                hold on
                plot(timeBar(:,1), timeBar(:,2), '--', 'linewidth', 1, 'color', style.fontColor);
                hold off
                title('Cell Voltage  /  V', 'color', style.fontColor)
                xlabel('Time  /  h')
                xlim([min(time/hour), max(time/hour)]);
                ylim([min(Enew), max(Enew)]);
                set(gca, ...
                    'FontSize' , style.fontSize       , ...
                    'FontName' , style.fontName       , ...
                    'color'    , style.backgroundColor, ...
                    'XColor'   , style.fontColor      , ...
                    'YColor'   , style.fontColor      , ...
                    'GridColor', style.fontColor)

            else
                style = setFigureStyle('theme', p.Results.theme, 'size', p.Results.size, 'orientation', p.Results.orientation, 'quantity', 'single');
                style.fontSize = 10;

                switch model.(ne).(am).diffusionModelType
                  case 'simple'
                    subplot(2,4,1), plotCellData(model.(ne).(am).G, states{i}.(ne).(am).c ./ 1000, 'edgealpha', 0.1);
                  case 'full'
                    subplot(2,4,1), plotCellData(model.(ne).(am).G, states{i}.(ne).(am).(sd).cSurface ./ 1000, 'edgealpha', 0.1);
                  otherwise
                    error('diffusionModelType not recognized');
                end
                xlabel(gca, 'Position  /  m')
                ylabel(gca, 'Position  /  m')
                title(gca, 'Negative Electrode Concentration  /  mol \cdot L^{-1}')
                colormap(crameri('nuuk'));
                caxis([cmin_ne, cmax_ne])
                colorbar
                axis tight
                set(gca, ...
                    'FontSize' , style.fontSize       , ...
                    'FontName' , style.fontName       , ...
                    'color'    , style.backgroundColor, ...
                    'XColor'   , style.fontColor      , ...
                    'YColor'   , style.fontColor      , ...
                    'GridColor', style.fontColor)
                if model.G.griddim == 3
                    view(45,45);
                    axis equal
                end


                subplot(2,4,2), plotCellData(model.(elyte).G, states{i}.(elyte).c ./ 1000, 'edgealpha', 0.1);
                xlabel(gca, 'Position  /  m')
                title(gca, 'Electrolyte Concentration  /  mol \cdot L^{-1}')
                colormap(crameri('nuuk'));
                caxis([cmin_elyte, cmax_elyte])
                colorbar
                axis tight
                set(gca, ...
                    'FontSize' , style.fontSize       , ...
                    'FontName' , style.fontName       , ...
                    'color'    , style.backgroundColor, ...
                    'XColor'   , style.fontColor      , ...
                    'YColor'   , style.fontColor      , ...
                    'GridColor', style.fontColor)
                if model.G.griddim == 3
                    view(45,45);
                    axis equal
                end

                switch model.(pe).(am).diffusionModelType
                  case 'simple'
                    subplot(2,4,3), plotCellData(model.(pe).(am).G, states{i}.(pe).(am).c ./ 1000, 'edgealpha', 0.1);
                  case 'full'
                    subplot(2,4,3), plotCellData(model.(pe).(am).G, states{i}.(pe).(am).(sd).cSurface ./ 1000, 'edgealpha', 0.1);
                  otherwise
                    error('diffusionModelType not recognized');
                end

                xlabel(gca, 'Position  /  m')
                title(gca, 'Positive Electrode Concentration  /  mol \cdot L^{-1}')
                colormap(crameri('nuuk'));
                caxis([cmin_pe, cmax_pe])
                colorbar
                axis tight
                set(gca, ...
                    'FontSize' , style.fontSize       , ...
                    'FontName' , style.fontName       , ...
                    'color'    , style.backgroundColor, ...
                    'XColor'   , style.fontColor      , ...
                    'YColor'   , style.fontColor      , ...
                    'GridColor', style.fontColor)
                if model.G.griddim == 3
                    view(45,45);
                    axis equal
                end

                subplot(2,4,4), plot((time/hour), Inew, '-', 'linewidth', 3)
                hold on
                plot(timeBar(:,1), timeBar(:,2), 'k--', 'linewidth', 1);
                hold off
                title('Cell Current  /  A')
                xlabel('Time  /  h')
                xlim([min(time/hour), max(time/hour)]);
                ylim([min(Inew), max(Inew)]);
                set(gca, ...
                    'FontSize' , style.fontSize       , ...
                    'FontName' , style.fontName       , ...
                    'color'    , style.backgroundColor, ...
                    'XColor'   , style.fontColor      , ...
                    'YColor'   , style.fontColor      , ...
                    'GridColor', style.fontColor)

                subplot(2,4,5), plotCellData(model.(ne).(am).G, states{i}.(ne).(am).phi, 'edgealpha', 0.1);
                xlabel(gca, 'Position  /  m')
                title(gca, 'Negative Electrode Potential  /  V')
                colormap(gca, crameri('lapaz'))
                caxis([phimin_ne, phimax_ne])
                colorbar
                axis tight
                set(gca, ...
                    'FontSize' , style.fontSize       , ...
                    'FontName' , style.fontName       , ...
                    'color'    , style.backgroundColor, ...
                    'XColor'   , style.fontColor      , ...
                    'YColor'   , style.fontColor      , ...
                    'GridColor', style.fontColor)
                if model.G.griddim == 3
                    view(45,45);
                    axis equal
                end

                subplot(2,4,6), plotCellData(model.(elyte).G, states{i}.(elyte).phi, 'edgealpha', 0.1);
                xlabel(gca, 'Position  /  m')
                title(gca, 'Electrolyte Potential  /  V')
                colormap(gca, crameri('lapaz'))
                caxis([phimin_elyte, phimax_elyte])
                colorbar
                axis tight
                set(gca, ...
                    'FontSize' , style.fontSize       , ...
                    'FontName' , style.fontName       , ...
                    'color'    , style.backgroundColor, ...
                    'XColor'   , style.fontColor      , ...
                    'YColor'   , style.fontColor      , ...
                    'GridColor', style.fontColor)
                if model.G.griddim == 3
                    view(45,45);
                    axis equal
                end

                subplot(2,4,7), plotCellData(model.(pe).(am).G, states{i}.(pe).(am).phi, 'edgealpha', 0.1);
                xlabel(gca, 'Position  /  m')
                title(gca, 'Positive Electrode Potential  /  V')
                colormap(gca, crameri('lapaz'))
                caxis([phimin_pe, phimax_pe])
                colorbar
                axis tight
                set(gca, ...
                    'FontSize' , style.fontSize       , ...
                    'FontName' , style.fontName       , ...
                    'color'    , style.backgroundColor, ...
                    'XColor'   , style.fontColor      , ...
                    'YColor'   , style.fontColor      , ...
                    'GridColor', style.fontColor)
                if model.G.griddim == 3
                    view(45,45);
                    axis equal
                end

                subplot(2,4,8), plot((time/hour), Enew, '-', 'linewidth', 3)
                hold on
                plot(timeBar(:,1), timeBar(:,2), 'k--', 'linewidth', 1);
                hold off
                title('Cell Voltage  /  V')
                xlabel('Time  /  h')
                xlim([min(time/hour), max(time/hour)]);
                ylim([min(Enew), max(Enew)]);
                set(gca, ...
                    'FontSize' , style.fontSize       , ...
                    'FontName' , style.fontName       , ...
                    'color'    , style.backgroundColor, ...
                    'XColor'   , style.fontColor      , ...
                    'YColor'   , style.fontColor      , ...
                    'GridColor', style.fontColor)
            end
            
            drawnow
            pause(0.1)
        end
    end
    
    
end

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
