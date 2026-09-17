classdef PEMgridGenerator1D < PEMgridGenerator


    properties

        % Length
        xlength
        % Discretization number
        N
        % Section area
        faceArea 
        
    end
    
    methods

        function gen = PEMgridGenerator1D()
            
            gen = gen@PEMgridGenerator();
            
        end

        function [inputparams, gen] = updateInputParams(gen, inputparams, params)

            inputparams = gen.setupInputParams(inputparams, []);
            
        end

        function [inputparams, gen] = setupGrid(gen, inputparams, params)
            
            G = cartGrid(gen.N, gen.xlength);

            parentGrid = Grid(G, 'faceArea', gen.faceArea);
            
            G = genSubGrid(parentGrid, (1 : parentGrid.getNumberOfCells())');
            
            inputparams.G = G;
            gen.parentGrid = parentGrid;
            
        end

        function inputparams = setupElectrodeElectrolyteCoupTerm(gen, inputparams, params)

            an    = 'Anode';
            ct    = 'Cathode';
            elyte = 'Electrolyte';

            G = gen.parentGrid;
            
            params.(an).couplingcells = 1;
            params.(an).couplingfaces = 1;
            params.(ct).couplingcells = G.getNumberOfCells();
            params.(ct).couplingfaces = G.getNumberOfFaces();

            inputparams = setupElectrodeElectrolyteCoupTerm@PEMgridGenerator(gen, inputparams, params);
            
        end

        function inputparams = setupControlArea(gen, inputparams)
        % Here inputparams is instance of ProtonicMembraneControlInputParams
            
            inputparams.area = gen.faceArea;
            
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
