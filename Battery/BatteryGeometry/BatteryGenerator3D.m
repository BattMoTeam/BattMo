classdef BatteryGenerator3D < BatteryGenerator
% Setup 3D grid with tab
    
    properties
            
        % Physical dimension
        
        xlength = 1e-2*[0.4; 0.2; 0.4];
        ylength = 1e-2*[0.1; 2; 0.1];
        zlength = 1e-6*[10; 100; 50; 80; 10];

        % Input current
        
        I = 1e-4;

        % Shorthands used below
        % ne    : NegativeElectrode
        % pe    : NegativeElectrode
        % am   : ElectrodActiveComponent
        % cc    : CurrentCollector
        % elyte : Electrolyte
        
        % Discretization resolution in z-direction
        
        facz = 1;
        
        sep_nz    = 3;
        ne_am_nz = 3;
        pe_am_nz = 3;
        ne_cc_nz  = 2;
        pe_cc_nz  = 2;
        
        % Discretization resolution in x-direction
        
        facx = 1;
        
        int_elyte_nx  = 3; 
        ne_cc_nx = 3;
        pe_cc_nx = 3;

        % Discretization resolution in y-direction

        facy = 1;
        
        ne_cc_ny = 2;
        pe_cc_ny = 2;
        elyte_ny = 4;

        % Utility variables computed once and then shared by methods (should not be set)
        elyte_nz;
        allparams;
        invcellmap;
        
        % Heat parameters
        externalHeatTransferCoefficientTab = 1e3;
        externalHeatTransferCoefficient = 1e3;

        use_thermal
                
    end
    
    methods
        
        function gen = BatteryGenerator3D()
            gen = gen@BatteryGenerator();  
        end
        
        function paramobj = updateBatteryInputParams(gen, paramobj)
            
            assert(paramobj.include_current_collectors, 'This geometry includes current collectors and input data appears to be missing for those');
            gen.use_thermal = paramobj.use_thermal;
            paramobj = gen.setupBatteryInputParams(paramobj, []);
            
        end
        
        function [paramobj, gen] = setupGrid(gen, paramobj, params)
            
            % shorthands
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            am    = 'ActiveMaterial';
            cc    = 'CurrentCollector';
            elyte = 'Electrolyte';
            sep   = 'Separator';
            
            gen = gen.applyResolutionFactors();
            
            nxs = [gen.ne_cc_nx; gen.int_elyte_nx; gen.pe_cc_nx];
            nys = [gen.ne_cc_ny; gen.elyte_ny; gen.pe_cc_ny];
            nzs = [gen.ne_cc_nz; gen.ne_am_nz; gen.sep_nz; gen.ne_am_nz; gen.pe_cc_nz];

            x = gen.xlength./nxs;
            x = rldecode(x, nxs);
            x = [0; cumsum(x)];

            y = gen.ylength./nys;
            y = rldecode(y, nys);
            y = [0; cumsum(y)];

            z = gen.zlength./nzs;
            z = rldecode(z, nzs);
            z = [0; cumsum(z)];

            G = tensorGrid(x, y, z);

            nx = sum(nxs);
            ny = sum(nys);
            nz = sum(nzs);

            dimGrid = [nx; ny; nz];
            
            gen.elyte_nz = gen.sep_nz + gen.ne_am_nz + gen.pe_am_nz;
            
            startSubGrid = [1; gen.ne_cc_ny + 1; gen.ne_cc_nz + 1];
            dimSubGrid   = [nx; gen.elyte_ny; gen.elyte_nz];
            allparams.(elyte).cellind = pickTensorCells3D(startSubGrid, dimSubGrid, dimGrid);
            
            startSubGrid = [1; gen.ne_cc_ny + 1; gen.ne_cc_nz + gen.ne_am_nz + 1];
            dimSubGrid   = [nx; gen.elyte_ny; gen.sep_nz];
            allparams.(elyte).(sep).cellind = pickTensorCells3D(startSubGrid, dimSubGrid, dimGrid);
            
            %% setup gen.ne_eac

            startSubGrid = [1; gen.ne_cc_ny + 1; gen.ne_cc_nz + 1];
            dimSubGrid   = [nx; gen.elyte_ny; gen.ne_am_nz];
            allparams.(ne).(am).cellind = pickTensorCells3D(startSubGrid, dimSubGrid, dimGrid);

            %% setup gen.pe_eac

            startSubGrid = [1; gen.ne_cc_ny + 1; gen.ne_cc_nz + gen.ne_am_nz + gen.sep_nz + 1];
            dimSubGrid   = [nx; gen.elyte_ny; gen.pe_am_nz];
            allparams.(pe).(am).cellind = pickTensorCells3D(startSubGrid, dimSubGrid, dimGrid);

            %% setup gen.ne_cc

            startSubGrid = [1; gen.ne_cc_ny + 1; 1];
            dimSubGrid   = [nx; gen.elyte_ny; gen.ne_cc_nz];
            allparams.(ne).(cc).cellind1 = pickTensorCells3D(startSubGrid, dimSubGrid, dimGrid);

            % We add the tab

            startSubGrid = [1; 1; 1];
            dimSubGrid   = [gen.ne_cc_nx; gen.ne_cc_ny; gen.ne_cc_nz];
            allparams.(ne).(cc).cellindtab = pickTensorCells3D(startSubGrid, dimSubGrid, dimGrid);

            allparams.(ne).(cc).cellind = [allparams.(ne).(cc).cellind1; allparams.(ne).(cc).cellindtab];

            %% setup gen.pe_cc

            startSubGrid = [1; gen.ne_cc_ny + 1; gen.ne_cc_nz + gen.elyte_nz + 1];
            dimSubGrid   = [nx; gen.elyte_ny; gen.pe_cc_nz];
            allparams.(pe).(cc).cellind1 = pickTensorCells3D(startSubGrid, dimSubGrid, dimGrid);

            % We add the tab

            startSubGrid = [gen.ne_cc_nx + gen.int_elyte_nx + 1; gen.ne_cc_ny + gen.elyte_ny + 1; gen.ne_cc_nz + gen.elyte_nz + 1];
            dimSubGrid   = [gen.pe_cc_nx; gen.pe_cc_ny; gen.pe_cc_nz];
            allparams.(pe).(cc).cellindtab = pickTensorCells3D(startSubGrid, dimSubGrid, dimGrid);

            allparams.(pe).(cc).cellind = [allparams.(pe).(cc).cellind1; allparams.(pe).(cc).cellindtab];

            cellind = [allparams.(elyte).cellind; allparams.(ne).(am).cellind; allparams.(pe).(am).cellind; allparams.(ne).(cc).cellind; allparams.(pe).(cc).cellind];

            rcellind = setdiff((1 : G.cells.num)', cellind);

            nGlob = G.cells.num;
            [G, cellmap, facemap, nodemap] = removeCells(G, rcellind);
            invcellmap = zeros(nGlob, 1);
            invcellmap(cellmap) = (1 : G.cells.num)';

            G = computeGeometry(G);
            
            paramobj.G = G;
            
            gen.invcellmap = invcellmap;
            gen.allparams = allparams;
            gen.G = G;

        end
        
        function gen = applyResolutionFactors(gen)
            
            facz = gen.facz;
            
            gen.sep_nz    = facz*gen.sep_nz;
            gen.ne_am_nz = facz*gen.ne_am_nz;
            gen.pe_am_nz = facz*gen.pe_am_nz;
            gen.ne_cc_nz  = facz*gen.ne_cc_nz;
            gen.pe_cc_nz  = facz*gen.pe_cc_nz;
            
            facx = gen.facx;
           
            gen.int_elyte_nx = facx*gen.int_elyte_nx;
            gen.ne_cc_nx= facx*gen.ne_cc_nx;
            gen.pe_cc_nx= facx*gen.pe_cc_nx;

            facy = gen.facy;
            
            gen.ne_cc_ny= facy*gen.ne_cc_ny;
            gen.pe_cc_ny= facy*gen.pe_cc_ny;
            gen.elyte_ny= facy*gen.elyte_ny;
            
        end

        function paramobj = setupElectrolyte(gen, paramobj, params)
            
            params = gen.allparams.Electrolyte;
            imap = gen.invcellmap;
            params.cellind = imap(params.cellind);
            params.Separator.cellind = imap(params.Separator.cellind);
            
            paramobj = setupElectrolyte@BatteryGenerator(gen, paramobj, params);
            
        end
        
        function paramobj = setupElectrodes(gen, paramobj, params)

            
            % shorthands 
            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            cc  = 'CurrentCollector';
            am = 'ActiveMaterial';

            params = gen.allparams;
            imap = gen.invcellmap;
            
            params.(ne).(am).cellind = imap(params.(ne).(am).cellind);
            params.(ne).(cc).cellind = imap(params.(ne).(cc).cellind);
            params.(ne).(cc).name = 'negative';
            params.(ne).cellind = [params.(ne).(am).cellind; params.(ne).(cc).cellind];
            
            params.(pe).(am).cellind = imap(params.(pe).(am).cellind);
            params.(pe).(cc).cellind = imap(params.(pe).(cc).cellind);
            params.(pe).(cc).name = 'positive';
            params.(pe).cellind = [params.(pe).(am).cellind; params.(pe).(cc).cellind];
            
            paramobj = setupElectrodes@BatteryGenerator(gen, paramobj, params);            
            
        end

        function paramobj = setupCurrentCollectorBcCoupTerm(gen, paramobj, params)
            
            G = paramobj.G;
            yf = G.faces.centroids(:, 2);
            
            switch params.name
              case 'negative'
                myf = min(yf);
              case 'positive'
                myf = max(yf);
            end
            
            params.bcfaces = find(abs(yf - myf) < eps*1000);
            params.bccells = sum(G.faces.neighbors(params.bcfaces, :), 2);

            paramobj = setupCurrentCollectorBcCoupTerm@BatteryGenerator(gen, paramobj, params);
        
        end

        function paramobj = setupThermalModel(gen, paramobj, params)
        % paramobj is instance of BatteryInputParams
        % 
        % We recover the external coupling terms for the current collectors

            % shorthands
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            cc    = 'CurrentCollector';
                
            % the cooling is done on the external faces
            G = gen.G;
            extfaces = any(G.faces.neighbors == 0, 2);
            couplingfaces = find(extfaces);
            couplingcells = sum(G.faces.neighbors(couplingfaces, :), 2);
            
            params = struct('couplingfaces', couplingfaces, ...
                            'couplingcells', couplingcells);
            paramobj = setupThermalModel@BatteryGenerator(gen, paramobj, params);
            
            tabcellinds = [gen.allparams.(pe).(cc).cellindtab; gen.allparams.(ne).(cc).cellindtab];
            tabtbl.cells = tabcellinds;
            tabtbl = IndexArray(tabtbl);
            
            tbls = setupSimpleTables(G);
            cellfacetbl = tbls.cellfacetbl;
            
            tabcellfacetbl = crossIndexArray(tabtbl, cellfacetbl, {'cells'});
            tabfacetbl = projIndexArray(tabcellfacetbl, {'faces'});
            
            bcfacetbl.faces = couplingfaces;
            bcfacetbl = IndexArray(bcfacetbl);
            
            tabbcfacetbl = crossIndexArray(bcfacetbl, tabfacetbl, {'faces'});
            
            map = TensorMap();
            map.fromTbl = bcfacetbl;
            map.toTbl = tabbcfacetbl;
            map.mergefds = {'faces'};
            ind = map.getDispatchInd();
            
            coef = gen.externalHeatTransferCoefficient*ones(bcfacetbl.num, 1);
            coef(ind) = gen.externalHeatTransferCoefficientTab;
            
            paramobj.ThermalModel.externalHeatTransferCoefficient = coef;
            
        end
        
    end
    
end


%{
Copyright 2021-2022 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
