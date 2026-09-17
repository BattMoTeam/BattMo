classdef FakeBatteryGeneratorP2D_doc < BatteryGeneratorP2D


    methods

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

            d = sqrt(gen.faceArea);
            y = [0; d];
            z = [0; d];
            
            G = tensorGrid(x, y, z);

            parentGrid = Grid(G, 'faceArea', gen.faceArea);

            G = genSubGrid(parentGrid, (1 : parentGrid.getNumberOfCells())');

            inputparams.G  = G;
            gen.parentGrid = parentGrid;


        end


    end

    
end

%{
Copyright 2021-2026 SINTEF Industry, Sustainable Energy Technology
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
