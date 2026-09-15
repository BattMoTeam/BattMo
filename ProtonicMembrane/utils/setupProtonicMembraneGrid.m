function [inputparams, gridGenerator] = setupProtonicMembraneGrid(inputparams, jsonstruct)

    an    = 'Anode';
    ct    = 'Cathode';
    elyte = 'Electrolyte';

    switch jsonstruct.Geometry.type

      case '1D'

        gen  = PEMgridGenerator1D();
        
        gen.xlength  = jsonstruct.(elyte).xlength;
        gen.faceArea = jsonstruct.(elyte).faceArea;
        gen.N        = jsonstruct.(elyte).Nx;

        inputparams = gen.updateInputParams(inputparams);

      case '2D'

        gen  = PEMgridGenerator2D();
        
        gen.ylength = jsonstruct.Geometry.ylength;
        gen.Ny      = jsonstruct.Geometry.Ny;
        gen.xlength = jsonstruct.(elyte).xlength;
        gen.Nx      = jsonstruct.(elyte).Nx;
        
        inputparams = gen.updateInputParams(inputparams);
        
      otherwise

        error('Geometry case not recognized');
        
    end
    
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
