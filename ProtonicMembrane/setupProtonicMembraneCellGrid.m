function paramobj = setupProtonicMembraneCellGrid(paramobj, jsonstruct)

    an    = 'Anode';
    ct    = 'Cathode';
    elyte = 'Electrolyte';

    lgth  = convertUnitBattMo(jsonstruct.(elyte).length);
    farea = jsonstruct.(elyte).faceArea;
    N     = jsonstruct.(elyte).N;

    dx = lgth/N;
    
    % setup grid

    G = cartGrid(N, lgth);
    G = computeGeometry(G);

    % Adjust for face area
    
    G.cells.volumes = farea*G.cells.volumes;
    G.faces.areas   = farea*G.faces.areas;
    G.faces.normals = farea*G.faces.normals;

    % Setup coupling terms

    couplingTerms = {};
    
    coupterm = couplingTerm('Anode-Electrolyte', {an, elyte});
    coupterm.couplingcells = [1, 1];
    coupterm.couplingfaces = [1, 1];
    couplingTerms{end + 1} = coupterm;
    
    coupterm = couplingTerm('Cathode-Electrolyte', {an, elyte});
    coupterm.couplingcells = [1, G.cells.num];
    coupterm.couplingfaces = [1, G.faces.num];
    couplingTerms{end + 1} = coupterm;

    % Assign to paramobj
    
    paramobj.G             = G;
    paramobj.farea         = farea;
    paramobj.dx            = dx;
    paramobj.(elyte).G     = G;
    paramobj.couplingTerms = couplingTerms;
    
end
