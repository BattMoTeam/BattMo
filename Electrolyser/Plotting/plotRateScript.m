setupoutput = false;
if setupoutput
    ind = cellfun(@(state) ~isempty(state), states);
    states = states(ind);
    for istate = 1 : numel(states)
        states{istate} = model.addVariables(states{istate});
    end
    time = cellfun(@(state) state.time, states);
    E = cellfun(@(state) state.(oer).(ptl).E, states);
    I = cellfun(@(state) state.(oer).(ctl).I, states);

    close all

    figure
    plot(time/hour, E)
    xlabel('time [hour]');
    ylabel('voltage');

    figure
    plot(time/hour, -I)
    xlabel('time [hour]');
    ylabel('Current [A]');
end

% load results from old code

close all

% str = 'H P ressure';
% str = ' P conc 2';
% str = 'H P j$';
str = 'P phi$';
% str = 'Cat Lay phi';
% str = 'phi$';
% str = 'cOHelyte';
% str = 'Cat Rate';
% str = 'Ox P B cOH';

% cgit.printVarNames(str);

% varinds = cgit.regexpVarNameSelect(str);

% fds = {'IonomerMembrane.j$'                                                        , ...
%        'OxygenEvolutionElectrode.PorousTransportLayer.H2OvaporLiquidExchangeRate'  , ...
%        'OxygenEvolutionElectrode.PorousTransportLayer.j$'                          , ...
%        'OxygenEvolutionElectrode.ExchangeReaction.H2OexchangeRate'                    , ...
%        'OxygenEvolutionElectrode.ExchangeReaction.OHexchangeRate'                     , ...
%        'OxygenEvolutionElectrode.CatalystLayer.inmrReactionRate$'                  , ...
%        'OxygenEvolutionElectrode.CatalystLayer.elyteReactionRate$'                 , ...
%        'HydrogenEvolutionElectrode.PorousTransportLayer.H2OvaporLiquidExchangeRate', ...
%        'HydrogenEvolutionElectrode.PorousTransportLayer.j$'                        , ...
%        'HydrogenEvolutionElectrode.ExchangeReaction.H2OexchangeRate'                  , ...
%        'HydrogenEvolutionElectrode.ExchangeReaction.OHexchangeRate'                   , ...
%        'HydrogenEvolutionElectrode.CatalystLayer.inmrReactionRate$'                , ...
%        'HydrogenEvolutionElectrode.CatalystLayer.elyteReactionRate$'};



fds = {'IonomerMembrane.j$'                                                                                                                       , ...
       'IonomerMembrane.OHsource'                                                                                                                 , ...
       'IonomerMembrane.H2OSource'                                                                                                                , ...
       'OxygenEvolutionElectrode.PorousTransportLayer.j$'                                                                                         , ...
       {'OxygenEvolutionElectrode.PorousTransportLayer.OHsource','OxygenEvolutionElectrode.PorousTransportLayer.OHbcSource'}                      , ...
       'OxygenEvolutionElectrode.PorousTransportLayer.H2OliquidSource'                                                                            , ...
       {'OxygenEvolutionElectrode.PorousTransportLayer.liquidSource','OxygenEvolutionElectrode.PorousTransportLayer.liquidBcSource'}              , ...
       {'OxygenEvolutionElectrode.PorousTransportLayer.compGasSources 1','OxygenEvolutionElectrode.PorousTransportLayer.compGasBcSources 1'}    , ...
       {'OxygenEvolutionElectrode.PorousTransportLayer.compGasSources 2','OxygenEvolutionElectrode.PorousTransportLayer.compGasBcSources 2'}    , ...
       'HydrogenEvolutionElectrode.PorousTransportLayer.j$'                                                                                       , ...
       {'HydrogenEvolutionElectrode.PorousTransportLayer.OHsource','HydrogenEvolutionElectrode.PorousTransportLayer.OHbcSource'}                  , ...
       'HydrogenEvolutionElectrode.PorousTransportLayer.H2OliquidSource'                                                                          , ...
       {'HydrogenEvolutionElectrode.PorousTransportLayer.liquidSource','HydrogenEvolutionElectrode.PorousTransportLayer.liquidBcSource'}          , ...
       {'HydrogenEvolutionElectrode.PorousTransportLayer.compGasSources 1','HydrogenEvolutionElectrode.PorousTransportLayer.compGasBcSources 1'}, ...
       {'HydrogenEvolutionElectrode.PorousTransportLayer.compGasSources 2','HydrogenEvolutionElectrode.PorousTransportLayer.compGasBcSources 2'}};

fdsdivs = {'IonomerMembrane.H2OdiffFlux'                                    , ...
           'IonomerMembrane.H2OmigFlux'                                     , ...
           'OxygenEvolutionElectrode.PorousTransportLayer.migOHFlux'        , ...
           'OxygenEvolutionElectrode.PorousTransportLayer.compGasFluxes 2'  , ...
           'OxygenEvolutionElectrode.PorousTransportLayer.compGasFluxes 1'  , ...
           'OxygenEvolutionElectrode.PorousTransportLayer.liquidMassFlux'   , ...
           'OxygenEvolutionElectrode.PorousTransportLayer.diffOHFlux'       , ...
           'OxygenEvolutionElectrode.PorousTransportLayer.convOHFlux'       , ...
           'HydrogenEvolutionElectrode.PorousTransportLayer.migOHFlux'      , ...
           'HydrogenEvolutionElectrode.PorousTransportLayer.compGasFluxes 2', ...
           'HydrogenEvolutionElectrode.PorousTransportLayer.compGasFluxes 1', ...
           'HydrogenEvolutionElectrode.PorousTransportLayer.liquidMassFlux' , ...
           'HydrogenEvolutionElectrode.PorousTransportLayer.diffOHFlux'     , ...
           'HydrogenEvolutionElectrode.PorousTransportLayer.convOHFlux'};


fds = horzcat(fds, fdsdivs);
% varinds = [];

% for ifd = 1 : numel(fds)
%     fd = fds{ifd};
%     ind = cgit.regexpVarNameSelect(fd);
%     assert(numel(ind) == 1, ['too many field selected for ' fd]);
%     varinds = [varinds; ind];
% end

% assert(numel(varinds) == numel(fds), 'setup problem : regexp was too greedy...');
dosave = true;

nt = numel(time);

for ivar = 1 : numel(fds)
    fd = fds{ivar};
    doaddbcterm = false;
    if iscell(fd)
        if doaddbcterm
            y = cell(nt, 1);
            for int = 1 : nt
                y{int} = 0;
            end
            for ifd = 1 : numel(fd)
                fdd = fd{ifd};
                varind = cgit.regexpVarNameSelect(fdd);
                assert(numel(varind) == 1, ['too many field selected for ' fdd]);
                varname = cgit.varNameList{varind};
                nodenames = cgit.getNodeName(varname);
                % should be unique here
                nnodename = nodenames{1};
                if ifd == 1
                    nodename = nnodename; % used later
                end
                [ymin, ymax, yy] = getvals(nnodename, states);
                for int = 1 : nt
                    y{int} = y{int} + yy{int};
                end
            end
        else
            fd = fd{1};
        end
    end

    if ~doaddbcterm
        varind = cgit.regexpVarNameSelect(fd);
        assert(numel(varind) == 1, ['too many field selected for ' fd]);
        varname = cgit.varNameList{varind};
        nodenames = cgit.getNodeName(varname);
        % should be unique here
        nodename = nodenames{1};
        [ymin, ymax, y] = getvals(nodename, states);
    end

    if ismember(fd, fdsdivs)
        dodiv = true;
    else
        dodiv = false;
    end

    if regexp(nodename, 'j$')
        dovol = false;
    else
        dovol = true;
    end

    if dodiv | dovol
        submodelname = regexp(nodename, '(.*)\.[^.]*$', 'tokens');
        submodelname = submodelname{1}{1};
        eval(sprintf('submodel = model.%s;', submodelname));
        vols = submodel.G.getVolumes();
    end

    figure
    hold on
    dominmax = false;
    if dominmax
        plot(time, ymin);
        plot(time, ymax);
    else
        imax = numel(time);
        for itime = 1 : imax
            if dodiv
                h = plot(submodel.G.getDiv(y{itime})./vols);
            elseif dovol
                h = plot(y{itime}./vols);
            else
                h = plot(y{itime});
            end
            if itime == 1
                set(h, 'linewidth', 3, 'color', 'b');
            elseif itime == imax
                set(h, 'linewidth', 3, 'color', 'r');
                if dodiv
                    title(sprintf('%s (div)', nodename), ' ');
                else
                    title(nodename, ' ');
                end
                if dosave
                    filename = sprintf('%s.png', nodename);
                    filename = fullfile(battmoDir, 'Electrolyser','img', filename);
                    saveas(gcf, filename);
                end
            end
        end
    end
end






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
