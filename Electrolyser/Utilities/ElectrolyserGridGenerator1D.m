classdef ElectrolyserGridGenerator1D < ElectrolyserGridGenerator


    properties

        xs % length of components in x direction
           % x(1) : length of oxygen porous transport layer
           % x(2) : length of oxygen catalyst layer
           % x(3) : length of pure ionomer (excluding the catalyst layer)
           % x(4) : length of hydrogen catalyst layer
           % x(5) : length of hydrogen porous transport layer

        cxs % length of the finite volume cells for each component (same dimension as )
        nxs % vector of same length as xs with the corresponding discretization numbers

    end

    methods

        function gen = ElectrolyserGridGenerator1D(xs, cxs)
            gen = gen@ElectrolyserGridGenerator();
            if nargin > 0
                gen.xs = xs;
                gen.cxs = cxs;
            end

        end

        function [inputparams, gen] = updateElectrolyserInputParams(gen, inputparams)

            xs  = gen.xs;
            cxs = gen.cxs;

            nxs = round(xs./cxs);
            gen.nxs = nxs;

            oer = 'OxygenEvolutionElectrode';
            her = 'HydrogenEvolutionElectrode';
            ptl = 'PorousTransportLayer';
            ctl = 'CatalystLayer';
            inm = 'IonomerMembrane';

            params.(inm).cellind = nxs(1) + (1 : sum(nxs(2 : 4)))';

            params.(oer).cellind       = (1 : sum(nxs(1 : 2)))';
            params.(oer).(ptl).cellind = (1 : sum(nxs(1 : 2)))';
            params.(oer).(ctl).cellind = nxs(1) + (1 : nxs(2))';
            params.(oer).bcfaces = 1;
            params.(oer).bccells = 1;

            params.(her).cellind       = sum(nxs(1 : 3)) + (1 : sum(nxs(4 : 5)))';
            params.(her).(ptl).cellind = sum(nxs(1 : 3)) + (1 : sum(nxs(4 : 5)))';
            params.(her).(ctl).cellind = sum(nxs(1 : 3)) + (1 : nxs(4))';
            params.(her).bcfaces = sum(nxs(4 : 5)) + 1; % index in own grid
            params.(her).bccells = sum(nxs(4 : 5)); % index in own grid

            [inputparams, gen] = gen.setupElectrolyserInputParams(inputparams, params);

        end


        function [inputparams, gen] = setupGrid(gen, inputparams, ~)

            xs = gen.xs;
            nxs = gen.nxs;

            x = xs./nxs;
            x = rldecode(x, nxs);
            x = [0; cumsum(x)];

            G = tensorGrid(x);

            parentGrid = Grid(G);

            G = genSubGrid(parentGrid, (1 : parentGrid.getNumberOfCells())');

            inputparams.G  = G;
            gen.parentGrid = parentGrid;

        end


    end


end



%{
Copyright 2021-2023 SINTEF Industry, Sustainable Energy Technology
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
