
doSetup = true;

if doSetup
    
    SOC = linspace(0, 1, 5);
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

%% Selected plots


newsimlist = bp.filterSimList(simlist, 'SOC', 1);
input = newsimlist{1};
output = outputs{input.simnumber};

%%

time = output.time;
E    = output.E;

figure
plot(time/year, E, '*-');
title('Voltage / V')
xlabel('Time / year')

%%

I = output.I;

figure
plot(time/year, I);
title('Current / A')
xlabel('Time / year')

%%

figure
hold on

plot(time/year, output.PE_Li_quantities                       ,'DisplayName', 'Positive Electrode');
plot(time/year, output.NE_Li_quantities                       ,'DisplayName', 'Negative Electrode');
plot(time/year, output.Electrolyte_Li_quantities              ,'DisplayName', 'Electrolyte');
plot(time/year, output.Electrodes_Li_quantities               ,'DisplayName', 'Both Electrodes');
plot(time/year, output.Total_Li_quantities                    ,'DisplayName', 'Total (except SEI)');
plot(time/year, output.Total_Li_quantities + output.quantities,'DisplayName', 'Total (including SEI)');
plot(time/year, output.quantities                             ,'DisplayName', 'In the SEI');

title('Lithium quantity');
xlabel('Time / year');
ylabel('quantity / mol');

grid on;

legend show

