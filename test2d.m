clear all
close all

%% add MRST module
mrstModule add ad-core

% setup object and run simulation
% delete('temp2.txt');
obj = lithiumIonModel();
obj.J = 1e-6;

%% run simulation

[t, y] = obj.p2d();

%% plotting
% Parse output

mrstModule add mrst-gui

fv = obj.fv;

compnames = obj.componentnames;

allstates = {};
Gs = {};

for icn = 1 : numel(compnames)
    
    compname = compnames{icn};
    comp = obj.(compname);

    Gs{icn} = comp.Grid;
    
    allstates{icn}= {};
    
    for iy = 1 : size(y, 1)
        yy = y(iy, :)';
        varname = 'phi';
        fullvarname = sprintf('%s_%s', compname, varname);
        state.(varname) = yy(fv.getSlot(fullvarname));
        if ismember(compname, {'ne', 'pe', 'elyte'})
            varname = 'Li';
            fullvarname = sprintf('%s_%s', compname, varname);
            state.(varname) = yy(fv.getSlot(fullvarname));
            if strcmp(compname, 'elyte')
               state.(varname) = state.(varname)./obj.elyte.eps;
            end
        end
        allstates{icn}{iy} = state;
    end
    
    figure(icn)
    plotToolbar(Gs{icn}, allstates{icn});
    title(compname);
    
end

% compnames = {'ne', 'ccne', 'pe', 'ccpe'};
% compnames = {'ne', 'ccne'};
compnames = {'pe', 'ccpe'};

states = {};
G = obj.G;
nc = G.cells.num;

for iy = 1 : size(y, 1)

    phi = nan(nc, 1);
    yy = y(iy, :)';
    
    for icn = 1 : numel(compnames)
    
        compname = compnames{icn};
        comp = obj.(compname);
        varname = 'phi';
        fullvarname = sprintf('%s_%s', compname, varname);
        philoc = yy(fv.getSlot(fullvarname));
        
        cellmap = comp.Grid.mappings.cellmap;
        
        phi(cellmap) = philoc;
        
    end

    states{iy}.phi = phi;
    
end

%%

figure
plotToolbar(G, states);
title('ne + ccne');


%% 


for iy = 1 : size(y, 1)
    yy = y(iy, :)';
    varname = 'E';
    E(iy) = yy(fv.getSlot(varname));
end

%%

figure
plot(E)





