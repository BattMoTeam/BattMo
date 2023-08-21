% load solution




elyte = 'Electrolyte';
ct    = 'Cathode';
ctam  = 'CathodeActiveMaterial';
an    = 'Anode';
anam  = 'AnodeActiveMaterial';

time = cellfun(@(x) x.time, states); 
Enew = cellfun(@(x) x.Cathode.E, states); 

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

dosavemovie = true;

if dosavemovie
    filename = 'volumefractions';
    avifilename = sprintf('%s.avi', filename);
    v = VideoWriter(avifilename);
    v.FrameRate = 3;
    open(v);
end


istate = 1;

doloop = true;
lastplottime = 0;
dplottime = 1*hour;
endtime = inf;

while doloop
    
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
    
    %% plot of species concentration
    % for ics  = 1 : 7
    %     if ~strcmp(model.(elyte).species(ics).name, 'H2O')
    %         plot(xelyte, log10(state.(elyte).cs{ics}/(mol/litre)), 'displayname', model.(elyte).species(ics).name);
    %     end
    % end
    % ics = 9;
    % plot(xelyte, log10(state.(elyte).cs{ics}/(mol/litre)), 'k-', 'displayname', model.(elyte).species(ics).name);
    
    %% plot of potential
    % val = state.Electrolyte.phi;
    % plot(xelyte, val, '*-', 'displayname', 'phi electrolyte');
    
    
    %% plot of quasi particles concentration
    % qpdict = model.(elyte).qpdict;
    % for key = qpdict.keys()
        % plot(xelyte, state.(elyte).qpcs{qpdict(key{1})}/(mol/litre), '*-', 'displayname', key{1});            
    % end
    
    % val = state.Electrolyte.Rprecipitation;
    % val = state.Electrolyte.volumeFraction;
    % val = state.Electrolyte.phi;
    % val = state.Electrolyte.eSource;
    % val = state.Electrolyte.solidVolumeFraction;
    % val = min(1, state.Electrolyte.nucleation);
    % val = state.Electrolyte.nucleation;
    % val = model.(elyte).G.getDiv(state.(elyte).qpFluxes{2});

    %% plot of precipitation rate
    % val = state.Electrolyte.Rprecipitation;
    % plot(xelyte, val, '*-', 'displayname', 'precipitation rate');
    
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
    
    legend('location', 'south west');
    xlabel('distance / m')
    if dosavemovie
        frame = getframe(gcf);
        writeVideo(v, frame);
    end
    drawnow
    pause(0.1)

end

if dosavemovie
    close(v)
    mp4filename = sprintf('%s.mp4', filename);
    cmd = sprintf('ffmpeg -y -i %s -vf \"crop=trunc(iw/2)*2:trunc(ih/2)*2\" %s', avifilename, mp4filename);
    system(cmd);
end

