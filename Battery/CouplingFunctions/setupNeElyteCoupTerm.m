function coupTerm = setupNeElyteCoupTerm(model)
    
    ne = model.ne;
    elyte = model.elyte;
    
    Gne = ne.G;
    Gelyte = elyte.G;
    
    % parent Grid
    G = Gne.mappings.parentGrid;
    
    % All the cells from ne are coupled with elyte
    cells1 = (1 : Gne.cells.num)';
    pcells = Gne.mappings.cellmap(cells1);
    
    mapping = zeros(G.cells.num, 1);
    mapping(Gelyte.mappings.cellmap) = (1 : Gelyte.cells.num)';
    cells2 = mapping(pcells);
    
    compnames = {'ne', 'elyte'};
    coupTerm = couplingTerm('ne-elyte', compnames);
    coupTerm.couplingcells =  [cells1, cells2];
    coupTerm.couplingfaces = []; % no coupling throug faces. We set it as empty
    
end
        