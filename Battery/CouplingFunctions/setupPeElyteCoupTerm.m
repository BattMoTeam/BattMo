function coupTerm = setupPeElyteCoupTerm(model)
    
    pe = model.pe;
    elyte = model.elyte;
    
    Gpe = pe.G;
    Gelyte = elyte.G;
    
    % parent Grid
    G = Gpe.mappings.parentGrid;
    
    % All the cells from pe are coupled with elyte
    cells1 = (1 : Gpe.cells.num)';
    pcells = Gpe.mappings.cellmap(cells1);
    
    mapping = zeros(G.cells.num, 1);
    mapping(Gelyte.mappings.cellmap) = (1 : Gelyte.cells.num)';
    cells2 = mapping(pcells);
    
    compnames = {'pe', 'elyte'};
    coupTerm = couplingTerm('pe-elyte', compnames);
    coupTerm.couplingcells = [cells1, cells2];
    coupTerm.couplingfaces = []; % no coupling between faces
    
end
