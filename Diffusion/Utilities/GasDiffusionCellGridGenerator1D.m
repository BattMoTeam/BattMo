classdef GasDiffusionCellGridGenerator1D < GasDiffusionCellGridGenerator


    properties

        length % length 
        N      % number of discretization cell

    end

    methods

        function gen = GasDiffusionCellGridGenerator1D(length, N)
            
            gen = gen@GasDiffusionCellGridGenerator();
            
            if nargin > 0
                gen.length = length;
                gen.N = N;
            end

        end

        function [inputparams, gen] = updateGasDiffusionCellInputParams(gen, inputparams, params)
        % this function is the main class function as it returns an updated inputparams object with grid structure

            if nargin < 3
                params = [];
            end
            [inputparams, gen] = setupGasDiffusionCellInputParams(gen, inputparams, params);
            
        end

        function [inputparams, gen] = setupGrid(gen, inputparams, ~)

            length = gen.length;
            N      = gen.N;

            G = cartGrid(N, length);

            parentGrid = Grid(G);

            G = genSubGrid(parentGrid, (1 : parentGrid.getNumberOfCells())');

            inputparams.G  = G;
            gen.parentGrid = parentGrid;

        end

        
        function inputparams = setupExternalCoupling(gen, inputparams, params)

            N = gen.N;

            params.nControls = 2;
            params.bcfaces = {[1], ...
                              [N + 1]};
            params.bccells = {[1], ...
                              [N]};

            inputparams = setupExternalCoupling@GasDiffusionCellGridGenerator(gen, inputparams, params);
            
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
