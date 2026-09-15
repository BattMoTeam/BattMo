function [inputparams, gridGenerator] = setupProtonicMembraneGasLayerGrid(inputparams, jsonstruct)

    geom = 'Geometry';
    
    switch jsonstruct.Geometry.type

      case '2D'

        gen = GasSupplyGridGenerator2D();

        gen.Nx      = jsonstruct.(geom).Nx;
        gen.Ny      = jsonstruct.(geom).Ny;
        gen.xlength = jsonstruct.(geom).xlength;
        gen.ylength = jsonstruct.(geom).ylength;
        
      case '1D'
        
        gen = GasSupplyGridGenerator1D();

        gen.Nx       = jsonstruct.(geom).Nx;
        gen.xlength  = jsonstruct.(geom).xlength;
        gen.faceArea = jsonstruct.(geom).faceArea;
        
      otherwise

        error('dimcase not recognized')
        
    end

    inputparams = gen.updateInputParams(inputparams);

    gridGenerator = gen;
    
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
