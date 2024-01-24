classdef SeaWaterBatteryGeneratorP2D < SeaWaterBatteryGenerator
% setup 1D grid
    properties

        andenx = 10; % number of cells in anode
        ctdenx = 10; % number of cells in cathode
        sepnx  = 10; % number of cells in electrolyte between the cathode and anode
        fac    = 1;

    end

    methods

        function gen = SeaWaterBatteryGeneratorP2D()
          gen = gen@SeaWaterBatteryGenerator();
        end

        function inputparams = updateBatteryInputParams(gen, inputparams)
            inputparams = gen.setupBatteryInputParams(inputparams, []);
        end

        function [inputparams, gen] = setupGrid(gen, inputparams, ~)
        % inputparams is instance of BatteryInputParams
        % setup inputparams.G

            andenx = gen.andenx;
            ctdenx = gen.ctdenx;
            sepnx  = gen.sepnx;

            nxs = [andenx; sepnx; ctdenx];

            % default value taken from original SeaWater code
            xlength = 1e-6*[5000; 1000; 100];

            x = xlength./nxs;
            x = rldecode(x, nxs);
            x = [0; cumsum(x)];

            G = tensorGrid(x);

            parentGrid = Grid(G);

            G = genSubGrid(parentGrid, (1 : parentGrid.getNumberOfCells())');

            inputparams.G  = G;
            gen.parentGrid = parentGrid;

        end

        function gen = applyResolutionFactors(gen)

            fac = gen.fac;

            gen.andenx = gen.andenx*fac;
            gen.ctdenx = gen.ctdenx*fac;
            gen.sepnx  = gen.sepnx*fac;

        end

        function inputparams = setupElectrolyte(gen, inputparams, params)

            % In this case we setup the electrolyte as a subgrid of the background, even if it the two in fact
            % coincides. It is unecessary but we do it to keep the approach generic.
            params.cellind = (1 : (gen.andenx + gen.sepnx + gen.ctdenx))';
            inputparams = setupElectrolyte@SeaWaterBatteryGenerator(gen, inputparams, params);
        end

        function inputparams = setupElectrodes(gen, inputparams, params)

            ande  = 'Anode';
            ctde  = 'Cathode';

            andenx = gen.andenx;
            ctdenx = gen.ctdenx;
            sepnx  = gen.sepnx;

            %% parameters for anode
            params.(ande).cellind = (1 : andenx)';
            params.(ande).bcfaces = 1;
            params.(ande).bccells = 1;

            %% parameters for cathode
            params.(ctde).cellind = (andenx + sepnx) + (1 : ctdenx)';
            params.(ctde).bcfaces = ctdenx + 1;
            params.(ctde).bccells = ctdenx;

            inputparams = setupElectrodes@SeaWaterBatteryGenerator(gen, inputparams, params);


        end            I

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
