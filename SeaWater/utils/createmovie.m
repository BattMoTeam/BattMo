function createmovie(casename, model, states, varargin)
    
    opt = struct('dosavemovie', true);
    opt = merge_options(opt, varargin{:});

    elyte = 'Electrolyte';
    ct    = 'Cathode';
    ctam  = 'CathodeActiveMaterial';
    an    = 'Anode';
    anam  = 'AnodeActiveMaterial';

    time = cellfun(@(x) x.time, states); 
    Enew = cellfun(@(x) x.Cathode.E, states); 

    
    if ismember('fullModel', fieldnames(model))
        model = model.fullModel;
    end
        
    
    ax = [min(time)/hour, max(time)/hour, min(Enew), max(Enew)];
    xelyte = model.Electrolyte.G.cells.centroids;
    xanode = model.Anode.G.cells.centroids;
    xcathode = model.Cathode.G.cells.centroids;
    qpdict = model.(elyte).qpdict;
    spdict = model.(elyte).spdict;

    set(0, 'defaulttextfontsize', 15);
    set(0, 'defaultaxesfontsize', 15);
    set(0, 'defaultlinelinewidth', 5);

    h = figure();
    set(h, 'position', [800, 600, 1600, 600]);
    figure(h)

    if opt.dosavemovie
        filename = casename;
        avifilename = sprintf('%s.avi', filename);
        v = VideoWriter(avifilename);
        v.FrameRate = 3;
        open(v);
    end

    istate = 1;

    doloop       = true;
    lastplottime = 0;
    dplottime    = 0;
    endtime      = inf;

    while doloop
        
        if dplottime == 0
            istate = istate + 1;
        end
        
            
        while (istate <= numel(states)) && (time(istate) - lastplottime < dplottime) 
            istate = istate + 1;
        end
        
        if istate > numel(states)
            istate = numel(states);
            doloop = false;
        end    

        if time(istate) > endtime
            doloop = false;
        end
        
        
        lastplottime = time(istate);
        
        subplot(1, 2, 1, 'replace')
        plot(time(1 : istate)/hour, Enew(1 : istate));
        xlabel('time / h')
        ylabel('potential at cathode / V');
        % title(sprintf('%g second, ind = %d', time(istate), istate));
        title(sprintf('%g second', time(istate)));
        axis(ax)

        subplot(1, 2, 2, 'replace')

        cla
        hold on
        state = states{istate};
        
        switch casename
          
          case 'concentrations'

            c = parula(floor(model.(elyte).nsp)/2);
            colororder(gca, c);
            set(gca, 'linestyleorder', {'-', '--'});
            
            %% plot of species concentration
            for ics  = 1 : model.(elyte).nsp
                if ~ismember(model.(elyte).species(ics).name, {'H2O', 'Na+'})
                    plot(xelyte, log10(state.(elyte).cs{ics}/(mol/litre)), 'displayname', model.(elyte).species(ics).name);
                end
            end
            
          case 'potential'

            %% plot of potential
            val = state.Electrolyte.phi;
            plot(xelyte, val, '*-', 'displayname', 'phi electrolyte');

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
        if opt.dosavemovie
            frame = getframe(gcf);
            writeVideo(v, frame);
        end
        drawnow
        pause(0.1)

    end

    if opt.dosavemovie
        close(v)
        mp4filename = sprintf('%s.mp4', filename);
        cmd = sprintf('ffmpeg -y -i %s -vf \"crop=trunc(iw/2)*2:trunc(ih/2)*2\" %s', avifilename, mp4filename);
        system(cmd);
    end

    
end
