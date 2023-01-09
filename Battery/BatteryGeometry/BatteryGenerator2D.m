classdef BatteryGenerator2D < BatteryGenerator
% Setup 2D grid

    properties

        xlength = 1e-6*[10; 100; 50; 80; 10]; % length of components in x direction - default values
                                              % x(1) : length of negative current collector
                                              % x(2) : length of negative active material
                                              % x(3) : length of separator
                                              % x(4) : length of positive active material
                                              % x(5) : length of positive current collector

        ylength = 1e-2; % length in y direction - default values

        ccnenx = 10; % discretization number for negative current collector - default value
        nenx   = 30; % discretization number for negative active material   - default value
        sepnx  = 30; % discretization number for separator                  - default value
        penx   = 30; % discretization number for positive active material   - default value
        ccpenx = 10; % discretization number for positive current collector - default value

        ny = 10; % discretization number in y direction - default value

        include_current_collectors % flag: true if grid for current collectors should be included

        use_thermal % flag, true if grid for thermal model should be setup.

        externalHeatTransferCoefficientTab = 1e3;  % heat transfer coefficient at tab boundary - default value
        externalHeatTransferCoefficient = 1e3;     % heat transfer coefficient at boundary - default value

    end

    methods

        function gen = BatteryGenerator2D()
            gen = gen@BatteryGenerator();
        end

        function paramobj = updateBatteryInputParams(gen, paramobj)

            fdnames = {'include_current_collectors', ...
                       'use_thermal'};

            gen = dispatchParams(gen, paramobj, fdnames);

            paramobj = gen.setupBatteryInputParams(paramobj, []);
        end

        function [paramobj, gen] = setupGrid(gen, paramobj, params)

            sepnx  = gen.sepnx;
            nenx   = gen.nenx;
            penx   = gen.penx;

            ylength = gen.ylength;

            if gen.include_current_collectors
                ccnenx = gen.ccnenx;
                ccpenx = gen.ccpenx;
                nxs = [ccnenx; nenx; sepnx; penx; ccpenx];
                xlength = gen.xlength;
            else
                nxs = [nenx; sepnx; penx];
                xlength = gen.xlength(2:4);
            end

            nx = sum(nxs);
            ny = gen.ny;

            x = xlength./nxs;
            x = rldecode(x, nxs);
            x = [0; cumsum(x)];

            y = ylength/ny;
            y = rldecode(y, ny);
            y = [0; cumsum(y)];

            G = tensorGrid(x, y);
            G = computeGeometry(G);

            paramobj.G = G;
            gen.G = G;

        end

        function paramobj = setupElectrolyte(gen, paramobj, params)

            sepnx  = gen.sepnx;
            nenx   = gen.nenx;
            penx   = gen.penx;

            if gen.include_current_collectors
                ccnenx = gen.ccnenx;
                ccpenx = gen.ccpenx;
                nxs = [ccnenx; nenx; sepnx; penx; ccpenx];
            else
                ccnenx = 0;
                ccpenx = 0;
                nxs = [nenx; sepnx; penx];
            end

            ny = gen.ny;
            nx = sum(nxs);

            istart = ccnenx + 1;
            ni = nenx + sepnx + penx;
            params.cellind = pickTensorCells(istart, ni, nx, ny);

            istart = ccnenx + nenx + 1;
            ni = sepnx;
            params.Separator.cellind = pickTensorCells(istart, ni, nx, ny);

            paramobj = setupElectrolyte@BatteryGenerator(gen, paramobj, params);

        end


        function paramobj = setupElectrodes(gen, paramobj, params)
        % setup grid and coupling term

            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            cc  = 'CurrentCollector';
            am = 'ActiveMaterial';

            sepnx  = gen.sepnx;
            nenx   = gen.nenx;
            penx   = gen.penx;

            if gen.include_current_collectors
                ccnenx = gen.ccnenx;
                ccpenx = gen.ccpenx;
                nxs = [ccnenx; nenx; sepnx; penx; ccpenx];
            else
                ccnenx = 0;
                ccpenx = 0;
                nxs = [nenx; sepnx; penx];
           end

            nx = sum(nxs);
            ny = gen.ny;

            %% Negative electrode
            istart = 1;
            ni = ccnenx + nenx;
            params.(ne).cellind = pickTensorCells(istart, ni, nx, ny);

            if gen.include_current_collectors
                % Negative electrode - current collector
                istart = 1;
                ni = ccnenx;
                params.(ne).(cc).cellind = pickTensorCells(istart, ni, nx, ny);
            end

            % Negative electrode - electrode active component
            istart = ccnenx + 1;
            ni = nenx;
            params.(ne).(am).cellind = pickTensorCells(istart, ni, nx, ny);

            %% Positive electrode
            istart = ccnenx + nenx + sepnx + 1;
            ni = penx + ccpenx;
            params.(pe).cellind = pickTensorCells(istart, ni, nx, ny);

            if gen.include_current_collectors
                % Positive electrode - current collector
                istart = ccnenx + nenx + sepnx + penx + 1;
                ni = ccpenx;
                params.(pe).(cc).cellind = pickTensorCells(istart, ni, nx, ny);
            end

            % Positive electrode - electrode active component
            istart = ccnenx + nenx + sepnx + 1;
            ni = penx;
            params.(pe).(am).cellind = pickTensorCells(istart, ni, nx, ny);

            if ~gen.include_current_collectors
                % Save electrode type for convenience
                params.(ne).(am).electrodeType = ne;
                params.(pe).(am).electrodeType = pe;
            end

            paramobj = setupElectrodes@BatteryGenerator(gen, paramobj, params);

        end

        function paramobj = setupCurrentCollectorBcCoupTerm(gen, paramobj, params)

            params = gen.findBoundary(paramobj.G, params);
            paramobj = setupCurrentCollectorBcCoupTerm@BatteryGenerator(gen, paramobj, params);

        end


        function paramobj = setupActiveMaterialBcCoupTerm(gen, paramobj, params)

            params = gen.findBoundary(paramobj.G, params);
            paramobj = setupActiveMaterialBcCoupTerm@BatteryGenerator(gen, paramobj, params);

        end


        function paramobj = setupThermalModel(gen, paramobj, params)
        % paramobj is instance of BatteryInputParams
        %
        % We recover the external coupling terms
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';

            if gen.include_current_collectors
                cc = 'CurrentCollector';
                coupterm_ne = paramobj.(ne).(cc).externalCouplingTerm;
                facemap_ne  = paramobj.(ne).(cc).G.mappings.facemap;
                coupterm_pe = paramobj.(pe).(cc).externalCouplingTerm;
                facemap_pe  = paramobj.(pe).(cc).G.mappings.facemap;
            else
                am = 'ActiveMaterial';
                coupterm_ne = paramobj.(ne).(am).externalCouplingTerm;
                facemap_ne  = paramobj.(ne).(am).G.mappings.facemap;
                coupterm_pe = paramobj.(pe).(am).externalCouplingTerm;
                facemap_pe  = paramobj.(pe).(am).G.mappings.facemap;
            end

            % the cooling is done on the external faces
            G = gen.G;
            extfaces = any(G.faces.neighbors == 0, 2);
            couplingfaces = find(extfaces);
            couplingcells = sum(G.faces.neighbors(couplingfaces, :), 2);

            params = struct('couplingfaces', couplingfaces, ...
                            'couplingcells', couplingcells);
            paramobj = setupThermalModel@BatteryGenerator(gen, paramobj, params);

            % We assign the values different values to the tab
            tabtbl.faces = [facemap_ne(coupterm_ne.couplingfaces);
                            facemap_pe(coupterm_pe.couplingfaces)];
            tabtbl = IndexArray(tabtbl);
            bcfacetbl.faces = couplingfaces;
            bcfacetbl = IndexArray(bcfacetbl);

            map = TensorMap();
            map.fromTbl = bcfacetbl;
            map.toTbl = tabtbl;
            map.mergefds = {'faces'};
            ind = map.getDispatchInd();

            coef = gen.externalHeatTransferCoefficient*ones(bcfacetbl.num, 1);
            coef(ind) = gen.externalHeatTransferCoefficientTab;

            paramobj.ThermalModel.externalHeatTransferCoefficient = coef;

        end

    end


    methods (Access = private)

        function params = findBoundary(gen, G, params)

            ne = 'NegativeElectrode';
            pe = 'PositiveElectrode';
            xf = G.faces.centroids(:, 1);

            switch params.electrodeType
              case ne
                x0 = min(xf);
              case pe
                x0 = max(xf);
            end

            params.bcfaces = find(abs(xf - x0) < eps*1000);
            params.bccells = sum(G.faces.neighbors(params.bcfaces, :), 2);

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
