function paramobj = setupBatteryInputParams1D(paramobj)
% params struct should contain fields
%
% - sepnx 
% - nenx 
% - penx 
% - ccnenx 
% - ccpenx 

    fac = 1;
    params.sepnx  = 30*fac;
    params.nenx   = 30*fac;
    params.penx   = 30*fac;
    params.ccnenx = 20*fac;
    params.ccpenx = 20*fac;
    
    paramobj = setupGrid(paramobj, params);
    paramobj = setupNegativeElectrode(paramobj, params);
    paramobj = setupPositiveElectrode(paramobj, params);
    paramobj = setupGridElectrolyte(paramobj, params);
    paramobj = setupCouplingTerms(paramobj, params);
    
end

function paramobj = setupGrid(paramobj, params)
    
    sepnx  = params.sepnx;
    nenx   = params.nenx;
    penx   = params.penx;
    ccnenx = params.ccnenx;
    ccpenx = params.ccpenx;

    nxs = [ccnenx; nenx; sepnx; penx; ccpenx];

    xlength = 1e-6*[10; 100; 50; 80; 10];
    ylength = 1e-2;

    x = xlength./nxs;
    x = rldecode(x, nxs);
    x = [0; cumsum(x)];

    G = tensorGrid(x);
    G = computeGeometry(G); 
    
end


function paramobj = setupNegativeElectrode(paramobj, params)
% setup grid and coupling term
    
    sepnx = params.sepnx; 
    nenx = params.nenx; 
    penx = params.penx; 
    ccnenx = params.ccnenx; 
    ccpenx = params.ccpenx;     
    
    params_ne.globG = paramobj.G;
    params_ne.cellind = (1 : ccnenx + nenx)';
    
    params_ne.eac.cellind = ccnenx + (1 : nenx)';
    params_ne.cc.cellind = (1 : ccnenx)';
    params_ne.cc.bc_cell = 1;
    params_ne.cc.bc_face = 1;
    
    params_ne.cc.eac_cell = ccnenx;
    params_ne.cc.eac_face = ccnenx + 1;
    params_ne.eac.cc_cell = 1;
    params_ne.eac.cc_face = 1;
    
    paramobj.ne = setupElectrodeInputParams1D(paramobj.ne, params_ne);
    
end

function pe_paramobj = setupPositiveElectrode(paramobj, params)
% setup grid and coupling term
    
    sepnx = params.sepnx; 
    nenx = params.nenx; 
    penx = params.penx; 
    ccnenx = params.ccnenx; 
    ccpenx = params.ccpenx;     
    
    pe_indstart = ccnenx + nenx + sepnx;
    
    params_pe.globG = paramobj.G;
    params_pe.cellind =  pe_indstart + (1 : ccpenx + penx)';
    
    params_pe.eac.cellind = pe_indstart + (1 : penx)';
    params_pe.cc.cellind = pe_indstart + penx + (1 : ccpenx)';
    params_pe.cc.bc_cell = ccpenx;
    params_pe.cc.bc_face = ccpenx + 1;
    params_pe.cc.eac_cell = 1;
    params_pe.cc.eac_face = 1;
    params_pe.eac.cc_cell = penx;
    params_pe.eac.cc_face = penx + 1;

    paramobj.pe = setupElectrodeInputParams1D(paramobj.pe, params_pe);
    
end

function elyte_paramobj = setupElectrolyte(paramobj, params)
% setup grid
    
    params_elyte.globG = paramobj.G;
    params_elyte.cellind = ccnenx + (1 : (nenx + sepnx + penx))';
    
    paramobj.elyte = setupGrid(paramobj.elyte, params_elyte);
   
end


function coupterms = setupCouplingTerms(paramobj, params)
    coupterms = {};
    coupterms{end + 1} = paramobj.setupCurrentCollectorBcCoupTerm(paramobj, param);
    coupterms{end + 1} = paramobj.setupCurrentCollectorElectrodeActiveComponentCoupTerm(paramobj, param);
end


function coupTerm = setupNegativeElectrodeElectrolyteCoupTerm(paramobj, params)
    
    nenx   = params.nenx;

    compnames = {'NegativeElectrode', 'Electrolyte'};
    coupTerm = couplingTerm('NegativeElectrode-Electrolyte', compnames);
    cells1 = (1 : nenx)';
    cells2 = (1 : nenx)';
    coupTerm.couplingcells =  [cells1, cells2];
    coupTerm.couplingfaces = []; % no coupling throug faces. We set it as empty
    
end

function coupTerm = setupPositiveElectrodeElectrolyteCoupTerm(paramobj, params)
    
    sepnx  = params.sepnx;
    nenx   = params.nenx;
    penx   = params.penx;
    
    compnames = {'PositiveElectrode', 'Electrolyte'};
    coupTerm = couplingTerm('PositiveElectrode-Electrolyte', compnames);
    cells1 = (1 : penx)';
    cells2 = nenx + sepnx + (1 : penx)';
    coupTerm.couplingcells = [cells1, cells2];
    coupTerm.couplingfaces = []; % no coupling between faces
    
end

