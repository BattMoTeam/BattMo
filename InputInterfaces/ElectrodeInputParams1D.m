classdef ElectrodeInputParams1D < ElectrodeInputParams
    
% Abbreviations used here:
% eac : ElectrodeActiveComponent
% cc  : CurrentCollector
% bc  : Boundary Condition
    
    methods
        
        function paramobj = ElectrodeInputParams1D(params)
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
        
            paramobj = paramobj@ElectrodeInputParams(params);
            
        end

        function eac_paramobj = setupElectrodeActiveComponent(paramobj, params)
            eac_param.globG = paramobj.G;
            eac_param.cellind = params.eac.cellind;
            
            eac_paramobj = ActiveElectroChemicalComponentInputParams(eac_params);
        end
        
        function cc_paramobj = setupCurrentCollector(paramobj, params)
            
            cc_params.globG = paramobj.G;
            cc_params.cellind = params.cc.cellind;
            cc_params.volumeFraction = params.cc.volumeFraction;
            cc_params.electronicConductivity = params.cc.electronicConductivity;
            
            cc_paramobj = CurrentCollectorInputParams(eac_params);
        
        end       
        
        function coupTerm = setupCurrentCollectorBcCoupTerm(paramobj, params)
        % Abbreviations used in this function:
        % cc   : CurrentCollector
            
            compnames = {'CurrentCollector'};
            coupTerm = couplingTerm('bc-CurrentCollector', compnames);
            coupTerm.couplingfaces = params.cc.bc_faces;
            coupTerm.couplingcells = params.cc.bc_cells;
            
        end
        
        function coupTerm = setupCurrentCollectorElectrodeActiveComponentCoupTerm(paramobj, params)
        % Abbreviations used in this function:
        % eac : ElectrodeActiveComponent
        % cc  : CurrentCollector
            
            compnames = {'CurrentCollector', 'ElectrodeActiveComponent'};
            coupTerm = couplingTerm('CurrentCollector-ElectrodeActiveComponent', compnames);
            coupTerm.couplingfaces =  [params.cc.eac_face, params.eac.cc_face];
            coupTerm.couplingcells =  [params.cc.eac_cell, params.eac.cc_cell];
            
        end

    end
    
end
