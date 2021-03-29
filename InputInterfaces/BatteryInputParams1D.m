classdef BatteryInputParams1D < BatteryInputParams
    
    methods
        
        function paramobj = BatteryInputParams1D(params)
        % params struct should contain valid fields for ComponentInputParams,
        %
        % and the fields
        %
        % - sepnx 
        % - nenx 
        % - penx 
        % - ccnenx 
        % - ccpenx 
        % - T
        % - SOC
        % - J
        % - Ucut
        % - ionName
        % - ionFluxName 
        % - ionSourceName
        % - ionMassConsName
        % - ionAccumName
            
            paramobj = paramobj@BatteryInputParams(params);
            
        end
        function G = setupGrid(paramobj, params)
            
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
        
        function ne_paramobj = setupNegativeElectrode(paramobj, params)
            
            sepnx = params.sepnx; 
            nenx = params.nenx; 
            penx = params.penx; 
            ccnenx = params.ccnenx; 
            ccpenx = params.ccpenx;     
            
            ne_params.globG = paramobj.G;
            ne_params.cellind = (1 : ccnenx + nenx)';
            
            ne_params.eac.cellind = ccnenx + (1 : nenx)';
            ne_params.cc.cellind = (1 : ccnenx)';
            ne_params.cc.bc_cell = 1;
            ne_params.cc.bc_face = 1;
            
            ne_params.cc.eac_cell = ccnenx;
            ne_params.cc.eac_face = ccnenx + 1;
            ne_params.eac.cc_cell = 1;
            ne_params.eac.cc_face = 1;

            ne_params.eac.ionName = params.ionName;
            ne_params.eac.ionFluxName  = params.ionFluxName ;
            ne_params.eac.ionSourceName = params.ionSourceName;
            ne_params.eac.ionMassConsName = params.ionMassConsName;
            ne_params.eac.ionAccumName = params.ionAccumName;
        
            ne_paramobj = ElectrodeInputParams1D(ne_params);
        
        end
        
        function pe_paramobj = setupPositiveElectrode(paramobj, params)
            
            sepnx = params.sepnx; 
            nenx = params.nenx; 
            penx = params.penx; 
            ccnenx = params.ccnenx; 
            ccpenx = params.ccpenx;     
            
            pe_indstart = ccnenx + nenx + sepnx;
            
            pe_params.globG = paramobj.G;
            pe_params.cellind =  pe_indstart + (1 : ccpenx + penx)';
            
            pe_params.eac.cellind = pe_indstart + (1 : penx)';
            pe_params.cc.cellind = pe_indstart + penx + (1 : ccpenx)';
            pe_params.cc.bc_cell = ccpenx;
            pe_params.cc.bc_face = ccpenx + 1;
            pe_params.cc.eac_cell = 1;
            pe_params.cc.eac_face = 1;
            pe_params.eac.cc_cell = penx;
            pe_params.eac.cc_face = penx + 1;

            pe_params.eac.ionName = params.ionName;
            pe_params.eac.ionFluxName  = params.ionFluxName ;
            pe_params.eac.ionSourceName = params.ionSourceName;
            pe_params.eac.ionMassConsName = params.ionMassConsName;
            pe_params.eac.ionAccumName = params.ionAccumName;
        
            pe_paramobj = ElectrodeInputParams1D(pe_params);

        end
        
        function elyte_paramobj = setupElectrolyte(paramobj, params)
            
            elyte_params.globG = paramobj.G;
            elyte_params.cellind = ccnenx + (1 : (nenx + sepnx + penx))';
            
            elyte_params.eac.ionName = params.ionName;
            elyte_params.eac.ionFluxName  = params.ionFluxName ;
            elyte_params.eac.ionSourceName = params.ionSourceName;
            elyte_params.eac.ionMassConsName = params.ionMassConsName;
            elyte_params.eac.ionAccumName = params.ionAccumName;
        
            elyte_paramobj = ElectrolyteInputParams(elyte_params);
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

    end
    
end
