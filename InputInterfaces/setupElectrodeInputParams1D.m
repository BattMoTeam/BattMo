function paramobj = setupElectrodeInputParams1D(paramobj, params)
    
% Abbreviations used here:
% eac : ElectrodeActiveComponent
% cc  : CurrentCollector
% bc  : Boundary Condition
    

% params struct should contain valid fields for ComponentInputParameters,
%
% and the fields
%
% - eac.cellind
% - cc.cellind
% - cc.bc_cell
% - cc.bc_face
% - cc.eac_cell
% - cc.eac_face
% - eac.cc_cell
% - eac.cc_face
    
    paramobj = paramobj.setupGrid(params);
    paramobj = setupElectrodeActiveComponent(paramobj, params);
    paramobj = setupCurrentCollector(paramobj, params);
    paramobj = setupCouplingTerms(paramobj, params);
    
end

function paramobj = setupElectrodeActiveComponent(paramobj, params)
    
    params_eac.globG = paramobj.globG;
    params_eac.cellind = params.eac.cellind;
    
    paramobj.eac = setupGrid(paramobj.eac, params_eac);
    
end

function paramobj = setupCurrentCollector(paramobj, params)
    
    params_cc.globG = paramobj.globG;
    params_cc.cellind = params.cc.cellind;
    
    paramobj.cc = setupGrid(paramobj.cc, params_cc);
    
end       

function paramobj = setupCouplingTerms(paramobj, params)
    
    coupterms = {};
    coupterms{end + 1} = setupCurrentCollectorBcCoupTerm(paramobj, params);
    coupterms{end + 1} = setupCurrentCollectorElectrodeActiveComponentCoupTerm(paramobj, params);
    paramobj.couplingTerms = coupterms;
    
end


function coupTerm = setupCurrentCollectorBcCoupTerm(paramobj, params)
% Abbreviations used in this function:
% cc   : CurrentCollector
    
    compnames = {'CurrentCollector'};
    coupTerm = couplingTerm('bc-CurrentCollector', compnames);
    coupTerm.couplingfaces = params.cc.bc_face;
    coupTerm.couplingcells = params.cc.bc_cell;
    
end

function coupTerm = setupCurrentCollectorElectrodeActiveComponentCoupTerm(paramobj, params)
% Abbreviations used in this function:
% eac : ElectrodeActiveComponent
% cc  : CurrentCollector
    
    compnames = {'CurrentCollector', 'ElectrodeActiveComponent'};
    coupTerm = couplingTerm('CurrentCollector-ElectrodeActiveComponent', compnames);
    coupTerm.couplingfaces = [params.cc.eac_face, params.eac.cc_face];
    coupTerm.couplingcells = [params.cc.eac_cell, params.eac.cc_cell];
    
end

