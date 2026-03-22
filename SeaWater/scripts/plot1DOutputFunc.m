function plot1DOutputFunc(simoutput, varargin)
%% plot 1D simulation output. The voltage curve until given time is plotted on a axes in the top left corner
%
% In varargin, casename gives the case to plot default is volumeFractions. The implemented "casename"s are
%   - 'concentrations'
%   - 'potential'
%   - 'quasiParticles'
%   - 'precipitationRate'
%   - 'nucleation'
%   - 'cSat'
%   - 'volumeFractions'
%   - 'elchemRates'
%
% In varargin, the figure number is also given
%
% The structure simoutput should have the fields that are obtained by calling getSimOutput and the ind for which the plot should be done, which are
%   - input  : input stucture as given as input for function run_magnesium_battery
%   - states : cell array of the output states
%   - model  : simulation model (MagnesiumBattery object as setup in run_magnesium_battery)
%   - ind    : index for states to plot
%
    
    opt = struct('casename', 'volumeFractions', ...
                 'fignum', 1);
    opt = merge_options(opt, varargin{:});
    
    elyte = 'Electrolyte';
    ct    = 'Cathode';
    ctam  = 'CathodeActiveMaterial';
    an    = 'Anode';
    anam  = 'AnodeActiveMaterial';

    fignum = opt.fignum;
    
    try
        set(0, 'currentFigure', fignum)
        allaxes = findall(gcf, 'type', 'axes');
        if ~(numel(allaxes) == 2 && strcmp(allaxes(1).Title.String, 'voltage'))
            clf
            error
        end            
    catch
        figure(fignum);
        set(fignum, 'Position', [100, 100, 2000, 1000]);
        axes('position', [0.1, 0.1, 0.8, 0.8]);
        axes('position', [0.05, 0.85, 0.1, 0.1]);
        allaxes = findall(gcf, 'type', 'axes');
    end

    input  = simoutput.input;
    ind    = simoutput.ind;
    states = simoutput.states;
    model  = simoutput.model;

    set(gcf, 'currentaxes', allaxes(1));
    
    time = cellfun(@(x) x.time, states); 
    Enew = cellfun(@(x) x.Cathode.E, states); 
    
    ax = [min(time)/hour, max(time)/hour, min(Enew), max(Enew)];
    
    plot(time(1 : ind)/hour, Enew(1 : ind));
    axis(ax)
    title('voltage')
    
    set(gcf, 'currentaxes', allaxes(2));
    
    cla
    hold on

    state = states{ind};
    
    xelyte   = model.Electrolyte.grid.cells.centroids(:, 1);
    xanode   = model.Anode.grid.cells.centroids(:, 1);
    xcathode = model.Cathode.grid.cells.centroids(:, 1);
    
    qpdict = model.(elyte).qpdict;
    spdict = model.(elyte).spdict;
    
    anvols    = model.(an).G.getVolumes();
    ctvols    = model.(an).G.getVolumes();
    elytevols = model.(elyte).G.getVolumes();
    F = model.con.F;
    
    switch opt.casename
        
      case 'concentrations'

        titlestr = 'Concentrations';
        c = parula(floor(model.(elyte).nsp));
        colororder(gca, c);
        set(gca, 'linestyleorder', {'-', '-.'});
        
        %% plot of species concentration
        for ics  = 1 : model.(elyte).nsp
            if ~ismember(model.(elyte).species(ics).name, {'H2O', 'Na+'})
                displayname = sprintf('%s',model.(elyte).species(ics).name);
                plot(xelyte, log10(state.(elyte).cs{ics}/(mol/litre)), 'displayname', displayname);
            end
        end
        
        ylabel('Concentrations / log10[mol / litre])');

      case 'potential'

        titlestr = 'Potential';
        
        %% plot of potential
        val = state.Electrolyte.phi;
        displayname =  'phi electrolyte';
        plot(xelyte, val, 'displayname', displayname);

        ylabel('Potential / [V]');

      case 'quasiParticles'

        titlestr = 'Quasi-particle concentrations';
        
        %% plot of quasi particles concentration
        vf = state.(elyte).volumeFraction;
        qpdict = model.(elyte).qpdict;
        for key = qpdict.keys()
            displayname = key{1};
            plot(xelyte, vf.*state.(elyte).qpepscs{qpdict(key{1})}/(mol/litre), 'displayname', displayname);            
        end
        
        ylabel('Concentrations / [mol / litre]');
        
      case 'precipitationRate'
        
        titlestr = 'Precipitation Rate';
        
        %% plot of precipitation rate
        val = state.Electrolyte.Rprecipitation;
        displayname =  'precipitation rate';
        plot(xelyte, val, 'displayname', displayname);
        
        ylabel('Rate / [mol m^-3 s^-1]');
        
      case 'nucleation'

        titlestr = 'Nucleation coefficient';
        
        %% plot of precipitation rate
        val = state.Electrolyte.nucleation;
        val = min(1, val);
        displayname =  'nucleation';
        plot(xelyte, val, 'displayname', displayname);
        
        ylabel('Nucleation coefficient / [-]');        
        
      case 'cSat'
        
        titlestr = 'Main ion concentration and saturation concentration';

        %% plot of precipitation rate
        c    = state.Electrolyte.cs{model.(elyte).mainIonIndex}/(mol/litre);
        csat = state.Electrolyte.cSat/(mol/litre);
        displayname =  'c';
        plot(xelyte, c, 'displayname', displayname);
        displayname =  'csat';
        plot(xelyte, csat, 'displayname', displayname);
        displayname =  'c - csat';
        plot(xelyte, c - csat, 'displayname', displayname);
        
        ylabel('[mol / litre]');


      case 'volumeFractions'
        
        titlestr = 'Volume fractions';

        %% plot of volume fractions
        val = state.Electrolyte.volumeFraction;
        displayname =  'electrolyte volume fraction';
        plot(xelyte, val, 'displayname', displayname);
        val = state.Electrolyte.solidVolumeFraction;
        displayname =  'solid volume fraction';
        plot(xelyte, val, 'displayname', displayname);
        val = state.(an).volumeFraction;
        displayname =  'Anode volume fraction';
        plot(xanode, val, 'displayname', displayname);
        val = state.(ct).volumeFraction;
        displayname = 'Cathode volume fraction';
        plot(xcathode, val, 'displayname', displayname);

              
        ylabel('Volume fractions / [-]');


      case 'elchemRates'
        
        titlestr = 'Redox reaction rates';

        %% plot of volume fractions
        val = state.(anam).R;
        displayname =  'Anode reaction rate';
        plot(xanode, val, 'displayname', displayname);
        val = state.(ctam).R;
        displayname = 'Cathode reaction rate';
        plot(xcathode, val, 'displayname', displayname);
              
        ylabel('Rates / [mol m^-3 s^-1])');
                  
      otherwise
        
        error('casename not recognized');
        
    end
    
    titlestr = sprintf('%s at time : %g hour', titlestr, state.time/hour);
    title(titlestr);
    legend('location', 'south west');
    xlabel('distance / m')

    drawnow

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
