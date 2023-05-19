function plotGivenTime(src, event, input)

    states    = input.states;
    model     = input.model;
    mainfigs  = input.mainfigs;
    voltagefig  = input.voltagefig;
    casenames = input.casenames;
    
    istate = get(src, 'Value');
    istate = floor(istate);
    
    elyte = 'Electrolyte';
    ct    = 'Cathode';
    ctam  = 'CathodeActiveMaterial';
    an    = 'Anode';
    anam  = 'AnodeActiveMaterial';

    time = cellfun(@(x) x.time, states); 
    Enew = cellfun(@(x) x.Cathode.E, states); 

    ax = [min(time)/hour, max(time)/hour, min(Enew), max(Enew)];

    set(0, 'currentFigure', voltagefig)
    
    plot(time(1 : istate)/hour, Enew(1 : istate));
    xlabel('time / h')
    ylabel('potential at cathode / V');
    % title(sprintf('%g second, ind = %d', time(istate), istate));
    title(sprintf('timestep %d, %g second', istate, time(istate)));
    axis(ax)    
    
    xelyte   = model.Electrolyte.G.cells.centroids;
    xanode   = model.Anode.G.cells.centroids;
    xcathode = model.Cathode.G.cells.centroids;
    qpdict   = model.(elyte).qpdict;
    spdict   = model.(elyte).spdict;
    
    for icase = 1 : numel(casenames)
        
        casename = casenames{icase};
        mainfig = mainfigs{icase};

        set(0, 'currentFigure', mainfig)
        
        cla
        hold on
        state = states{istate};
        
        switch casename

          case 'an_specific_area'

            %% plot of potential
            val = state.(anam).specificSurfaceArea;
            plot(xanode, val, '*-', 'displayname', 'anode specific area');
            
            
          case 'concentrations'

            c = parula(floor(model.(elyte).nsp)/2);
            colororder(gca, c);
            set(gca, 'linestyleorder', {'-', '-.'});
            
            %% plot of species concentration
            for ics  = 1 : model.(elyte).nsp
                % if ~ismember(model.(elyte).species(ics).name, {'H2O', 'Na+', 'Mg(OH)2'})
                if ~ismember(model.(elyte).species(ics).name, {'H2O', 'Na+'})
                    plot(xelyte, log10(state.(elyte).cs{ics}/(mol/litre)), 'displayname', model.(elyte).species(ics).name);
                end
            end
            
          case 'potential'

            %% plot of potential
            % val = state.Electrolyte.phi;
            % plot(xelyte, val, '*-', 'displayname', 'phi electrolyte');
            % val = state.(ctam).eta;
            % plot(xcathode, val, '*-', 'displayname', 'eta cathode');
            val = state.(anam).eta;
            plot(xanode, val, '*-', 'displayname', 'eta anode');        
            
          case 'precipitationSurfaceArea'

            %% plot of potential
            val = state.Electrolyte.precipitationSurfaceArea;
            plot(xelyte, val, '*-', 'displayname', 'precipitationSurfaceArea');
            
          case 'covercsat'
            
            %% plot of precipitation rate
            c    = state.Electrolyte.cs{model.(elyte).mainIonIndex};
            csat = state.Electrolyte.cSat;
            % val = c./csat;
            plot(xelyte, c, '*-', 'displayname', 'c');
            plot(xelyte, csat, '*-', 'displayname', 'csat');
            plot(xelyte, c- csat, '*-', 'displayname', 'c - csat');
            
          case 'quasiparticles'
            
            %% plot of quasi particles concentration
            qpdict = model.(elyte).qpdict;
            for key = qpdict.keys()
                plot(xelyte, state.(elyte).qpcs{qpdict(key{1})}/(mol/litre), '*-', 'displayname', key{1});            
            end
            
          case 'precipitationrate'
            
            %% plot of precipitation rate
            val = state.Electrolyte.Rprecipitation;
            plot(xelyte, val, '*-', 'displayname', 'precipitation rate');
            
          case 'nucleation'
            
            %% plot of precipitation rate
            val = min(1, state.Electrolyte.nucleation);
            plot(xelyte, val, '*-', 'displayname', 'nucleation');
            axis([xelyte(1), xelyte(end), -0.1, 1.5])
            
          case 'volumefractions'
            
            %% plot of volume fractions
            val = state.Electrolyte.volumeFraction;
            plot(xelyte, val, '*-', 'displayname', 'electrolyte volume fraction');
            val = state.Electrolyte.solidVolumeFraction;    
            plot(xelyte, val, '*-', 'displayname', 'solid volume fraction');
            val = state.(an).volumeFraction;
            plot(xanode, val, '*-', 'displayname', 'Anode volume fraction');
            val = state.(ct).volumeFraction;
            plot(xcathode, val, '*-', 'displayname', 'Cathode volume fraction');
            
            %% plot of solid concentration
            % val = log10(state.Electrolyte.cs{spdict('Mg(OH)2')}/(mol/litre));
            % plot(xelyte, val, '*-', 'displayname', 'log10(MgOH2)');

            
            %% plot of saturation
            % val = state.Electrolyte.nucleation;
            % plot(xelyte, min(1, val), '*-', 'displayname', 'nucleation flag');
            % val = log10(state.Electrolyte.cs{spdict('Mg+2')}/(mol/litre));
            % plot(xelyte, val, '*-', 'displayname', 'log(Mg+2)');    
            % val = log10(state.Electrolyte.cSat/(mol/litre));
            % plot(xelyte, val, '*-', 'displayname', 'log((Mg+2)_{sat})');    
            
            %% plot of values at anode
            % val = state.(anam).ENernst;
            % plot(xanode, val, '*-', 'displayname', 'Enernst potential - Anode');
            % val = state.(ctam).ENernst;
            % plot(xcathode, val, '*-', 'displayname', 'Enernst potiential - Cathode');
            % val = state.(an).phi;
            % plot(xanode, val, '*-', 'displayname', 'Anode phi');
            % val = state.(ct).phi;
            % plot(xcathode, val, '*-', 'displayname', 'Cathode  phi');
            % val = state.(anam).eta;
            % plot(xanode, val, '*-', 'displayname', 'Anode eta');
            % val = state.(ctam).eta;
            % plot(xcathode, val, '*-', 'displayname', 'Cathode  eta');
            % val = state.(elyte).phi;
            % plot(xelyte, val, '*-', 'displayname', 'Electrolyte phi');
            
            %% plot of values at cathode
            % val = log10(state.CathodeActiveMaterial.cElectrolyte);
            % plot(xcathode, val, '*-');
            
          otherwise
            
            error('casename not recognized');
            
        end
        
        legend('location', 'south west');
        xlabel('distance / m')
    end
end
