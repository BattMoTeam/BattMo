classdef BatteryGeneratorP2D < BatteryGenerator
% setup 1D grid
    properties

        %
        % vector of component lengths
        %
        % - x(1) : length of negative current collector (default = 25 micro meter)
        % - x(2) : length of negative active material (default = 64 micro meter)
        % - x(3) : length of separator (default = 15 micro meter)
        % - x(4) : length of positive active material (default = 57 micro meter)
        % - x(5) : length of positive current collector (default = 15 micro meter)
        %
        xlength = 1e-6*[25; 64; 15; 57; 15];

        sepnx  = 10; % discretization number for negative current collector (default = 10)
        nenx   = 10; % discretization number for negative active material (default = 10)
        penx   = 10; % discretization number for separator (default = 10)
        ccnenx = 10; % discretization number for positive current collector (default = 10)
        ccpenx = 10; % discretization number for positive active material (default = 10)

        %
        % refinement factor (can be used to easily increase discretization refinement)
        % see applyResolutionFactors method
        %
        resolutionFactor = 1;

        % boolean : true if grid for current collectors should be included
        include_current_collectors
        % boolean : true if grid for thermal model should be included
        use_thermal

        % Face area in the transversal direction (default = 1)
        faceArea = 2*1.6387e-04;

    end

    methods


        function gen = BatteryGeneratorP2D()
            gen = gen@BatteryGenerator();
        end


        function [inputparams, gen] = updateBatteryInputParams(gen, inputparams)

            gen.include_current_collectors = inputparams.include_current_collectors;
            gen.use_thermal = inputparams.use_thermal;
            inputparams = gen.setupBatteryInputParams(inputparams, []);

        end


        function [inputparams, gen] = setupGrid(gen, inputparams, ~)

            sepnx  = gen.sepnx;
            nenx   = gen.nenx;
            penx   = gen.penx;
            ccnenx = gen.ccnenx;
            ccpenx = gen.ccpenx;

            xlength = gen.xlength;

            if gen.include_current_collectors
                nxs = [ccnenx; nenx; sepnx; penx; ccpenx];
            else
                nxs = [nenx; sepnx; penx];
                xlength = xlength(2 : 4);
            end

            x = xlength./nxs;
            x = rldecode(x, nxs);
            x = [0; cumsum(x)];

            G = tensorGrid(x);

            parentGrid = Grid(G, 'faceArea', gen.faceArea);

            G = genSubGrid(parentGrid, (1 : parentGrid.getNumberOfCells())');

            inputparams.G  = G;
            gen.parentGrid = parentGrid;

        end


        function gen = applyResolutionFactors(gen)

            fac = gen.resolutionFactor;

            gen.sepnx  = gen.sepnx*fac;
            gen.nenx   = gen.nenx*fac;
            gen.penx   = gen.penx*fac;
            gen.ccnenx = gen.ccnenx*fac;
            gen.ccpenx = gen.ccpenx*fac;

        end


        function inputparams = setupElectrolyte(gen, inputparams, params)

            if gen.include_current_collectors
                params.cellind = gen.ccnenx + (1 : (gen.nenx + gen.sepnx + gen.penx))';
            else
                params.cellind = (1 : (gen.nenx + gen.sepnx + gen.penx))';
            end
            inputparams = setupElectrolyte@BatteryGenerator(gen, inputparams, params);

        end


        function inputparams = setupSeparator(gen, inputparams, params)

            if gen.include_current_collectors
                params.cellind = gen.ccnenx + gen.nenx + (1 : gen.sepnx)';
            else
                params.cellind =  gen.nenx + (1 : gen.sepnx)';
            end
            inputparams = setupSeparator@BatteryGenerator(gen, inputparams, params);

        end


        function inputparams = setupElectrodes(gen, inputparams, params)

            ne = 'NegativeElectrode';
            pe = 'PositiveElectrode';
            cc = 'CurrentCollector';
            co = 'Coating';

            sepnx  = gen.sepnx;
            nenx   = gen.nenx;
            penx   = gen.penx;

            if gen.include_current_collectors

                ccnenx = gen.ccnenx;
                ccpenx = gen.ccpenx;

                %% parameters for negative electrode

                params.(ne).cellind = (1 : ccnenx + nenx)';
                params.(ne).(co).cellind = ccnenx + (1 : nenx)';
                params.(ne).(cc).cellind = (1 : ccnenx)';

                % boundary setup for negative current collector
                params.(ne).(cc).bcfaces = 1;
                params.(ne).(cc).bccells = 1;

                %% parameters for positive electrode

                pe_indstart = ccnenx + nenx + sepnx;
                params.(pe).cellind =  pe_indstart + (1 : ccpenx + penx)';
                params.(pe).(co).cellind = pe_indstart + (1 : penx)';
                params.(pe).(cc).cellind = pe_indstart + penx + (1 : ccpenx)';

                % boundary setup for positive current collector
                params.(pe).(cc).bcfaces = ccpenx + 1;
                params.(pe).(cc).bccells = ccpenx;

            else

                %% parameters for negative electrode

                params.(ne).cellind = (1 : nenx)';
                params.(ne).(co).cellind = (1 : nenx)';

                % boundary setup for negative current collector
                params.(ne).(co).bcfaces = 1;
                params.(ne).(co).bccells = 1;

                %% parameters for positive electrode

                pe_indstart = nenx + sepnx;
                params.(pe).cellind =  pe_indstart + (1 :  penx)';
                params.(pe).(co).cellind = pe_indstart + (1 : penx)';

                % boundary setup for positive current collector
                params.(pe).(co).bcfaces = penx + 1;
                params.(pe).(co).bccells = penx;

            end

            inputparams = setupElectrodes@BatteryGenerator(gen, inputparams, params);

        end


        function inputparams = setupThermalModel(gen, inputparams, params)

            params.couplingfaces = [];
            params.couplingcells = (1 : gen.parentGrid.getNumberOfCells())';
            inputparams = setupThermalModel@BatteryGenerator(gen, inputparams, params);

        end

    end


end





%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
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
