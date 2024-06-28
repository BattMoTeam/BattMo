
doSetup = false;

if doSetup
    
    SOC = linspace(0, 1, 15);
    % DiffE = [2; 3; 4]*1e-15;
    % [~, simlist] = cartesianProduct(SOC, DiffE);
    
    [~, simlist] = cartesianProduct(SOC);
    
    outputs = {};
    for isim = 1 : numel(simlist)
        
        input = simlist{isim};
        simlist{isim}.simnumber = isim;
        outputs{isim} = runBolaySEIfunction(input);
        
    end

end

bp = BatchProcessor(simlist);
bp.printSimList(simlist);

capacity = computeCellCapacity(outputs{1}.model);
F = PhysicalConstants.F;

figure
hold on

for isim = 1 : numel(simlist)

    input = simlist{isim};

    output = outputs{input.simnumber};
    
    qty  = output.quantities;
    qty  = qty - qty(1);
    time = output.time;
    
    % legendtxt = sprintf('SOC : %g, De : %g', input.SOC, input.DiffE);
    legendtxt = sprintf('SOC : %g', input.SOC);
    plot(time/year, 100 * (qty*F) / capacity, 'displayname', legendtxt);
    
end

title('Percentage of Lithium consummed');
xlabel('Time / year');
ylabel('%');

legend show

grid on;

qtys = [];
socs = [];
for isim = 1 : numel(simlist)

    input = simlist{isim};

    socs(isim) = input.SOC;
    
    output = outputs{input.simnumber};
    
    qty = output.quantities;
    qty = qty - qty(1);
    qty = 100 * (qty*F) / capacity;
    
    qtys(isim) = qty(end);
    
end

figure
plot(socs, qtys)

return


%% Setup for plotting

ind = cellfun(@(x) not(isempty(x)), states); 
states = states(ind);
time = cellfun(@(x) x.time, states); 
E    = cellfun(@(x) x.Control.E, states); 
I    = cellfun(@(x) x.Control.I, states);

for istate = 1 : numel(states)
    states{istate} = model.addVariables(states{istate});
end

%% Plotting

close all

set(0, 'defaultlinelinewidth', 3)
set(0, 'defaultaxesfontsize', 15)

%%

figure
plot(time/year, E, '*-');
title('Voltage / V')
xlabel('Time / year')

%%

figure
plot(time/year, I);
title('Current / A')
xlabel('Time / year')

%%

figure
hold on

delta = cellfun(@(state) state.(ne).(co).(am).(itf).SEIlength(end), states);
plot(time/year, delta/(nano*meter), 'displayname', 'at x_{max}')

delta = cellfun(@(state) state.(ne).(co).(am).(itf).SEIlength(1), states);
plot(time/year, delta/(nano*meter), 'displayname', 'at x_{min}')

title('SEI thickness in negative electrode/ nm')
xlabel('Time / year')

legend show

%%

figure
hold on

u = cellfun(@(state) state.(ne).(co).(am).(itf).SEIvoltageDrop(end), states);
plot(time/year, u, 'displayname', 'at x_{max}')

u = cellfun(@(state) state.(ne).(co).(am).(itf).SEIvoltageDrop(1), states);
plot(time/year, u, 'displayname', 'at x_{min}')

title('SEI voltage drop in negative electrode/ V')
xlabel('Time / year')

%%
figure

vols = model.(ne).(co).G.getVolumes();

for istate = 1 : numel(states)
    state = states{istate};
    m(istate) = sum(vols.*state.(ne).(co).(am).(sd).cAverage);
end

plot(time/year, m);
title('total lithium amount in negative electrode / mol')
xlabel('Time / year')

legend show

%%

quantities = [];

vols = model.(ne).(co).G.getVolumes();

for timeindex = 1 : numel(states)

    state = states{timeindex};
    state = model.evalVarName(state, {ne, co, am, itf, 'SEIconcentration'});
    cSEI  = state.(ne).(co).(am).(itf).SEIconcentration;
    
    Liqqt = sum(cSEI.*vols);
    quantities(end + 1) = Liqqt;
        
end

figure 
plot(time/year, quantities);
title('Lithium quantity consummed');
xlabel('Time / year');
ylabel('quantity / mol');
grid on;

%%

PE_Li_quantities          = [];
NE_Li_quantities          = [];
Electrolyte_Li_quantities = [];
Electrodes_Li_quantities  = [];
Total_Li_quantities       = [];

for timeindex = 1 : numel(states)


    amvf     = model.(pe).(co).volumeFractions(1);
    vf       = model.(pe).(co).volumeFraction;
    vols     = model.(pe).G.getVolumes;
    cAverage = states{timeindex}.(pe).(co).(am).(sd).cAverage;

    PE_qtt = sum(amvf.*vf.*vols.*cAverage);
    
    amvf     = model.(ne).(co).volumeFractions(1);
    vf       = model.(ne).(co).volumeFraction;
    vols     = model.(ne).G.getVolumes;
    cAverage = states{timeindex}.(ne).(co).(am).(sd).cAverage;

    NE_qtt = sum(amvf.*vf.*vols.*cAverage);

    Elyte_qtt = sum(model.Electrolyte.volumeFraction.*model.Electrolyte.G.getVolumes.*states{timeindex}.Electrolyte.c);

    Elode_qtt = PE_qtt + NE_qtt;
    Tot_Liqqt = PE_qtt + NE_qtt + Elyte_qtt;
	
    PE_Li_quantities(end + 1)          = PE_qtt;
    NE_Li_quantities(end + 1)          = NE_qtt;
    Electrolyte_Li_quantities(end + 1) = Elyte_qtt;
    Electrodes_Li_quantities(end + 1)  = Elode_qtt;
    Total_Li_quantities(end + 1)       = Tot_Liqqt;

end

%%

figure
hold on

plot(time/year, PE_Li_quantities                ,'DisplayName', 'Positive Electrode');
plot(time/year, NE_Li_quantities                ,'DisplayName', 'Negative Electrode');
plot(time/year, Electrolyte_Li_quantities       ,'DisplayName', 'Electrolyte');
plot(time/year, Electrodes_Li_quantities        ,'DisplayName', 'Both Electrodes');
plot(time/year, Total_Li_quantities             ,'DisplayName', 'Total (except SEI)');
plot(time/year, Total_Li_quantities + quantities,'DisplayName', 'Total (including SEI)');
plot(time/year, quantities                      ,'DisplayName', 'In the SEI');

title('Lithium quantity');
xlabel('Time / year');
ylabel('quantity / mol');

grid on;

legend show

%%

capacity = computeCellCapacity(model);
F = PhysicalConstants.F;

figure 
plot(time/year, 100 * (quantities*F) / capacity);
title('Percentage of Lithium consummed');
xlabel('Time / year');
ylabel('%');
grid on;

