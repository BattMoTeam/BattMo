function cells = pickTensorCells3D(startSubGrid, dimSubGrid, dimGlobalGrid)
    
    Gxtbl.x = (1 : dimGlobalGrid(1))';
    Gytbl.y = (1 : dimGlobalGrid(2))';
    Gztbl.z = (1 : dimGlobalGrid(3))';

    Gxtbl = IndexArray(Gxtbl);
    Gytbl = IndexArray(Gytbl);
    Gztbl = IndexArray(Gztbl);    
    
    Gtbl = crossIndexArray(Gxtbl, Gytbl, {}, 'optpureproduct', true);
    Gtbl = crossIndexArray(Gtbl, Gztbl, {}, 'optpureproduct', true);
    
    Gtbl = sortIndexArray(Gtbl, {'z', 'y', 'x'});
    
    Lxtbl.x = (startSubGrid(1) : (startSubGrid(1) + dimSubGrid(1) - 1))';
    Lytbl.y = (startSubGrid(2) : (startSubGrid(2) + dimSubGrid(2) - 1))';
    Lztbl.z = (startSubGrid(3) : (startSubGrid(3) + dimSubGrid(3) - 1))';

    Lxtbl = IndexArray(Lxtbl);
    Lytbl = IndexArray(Lytbl);
    Lztbl = IndexArray(Lztbl);        
    
    Ltbl = crossIndexArray(Lxtbl, Lytbl, {}, 'optpureproduct', true);
    Ltbl = crossIndexArray(Ltbl, Lztbl, {}, 'optpureproduct', true);    
    
    Ltbl = sortIndexArray(Ltbl, {'z', 'y', 'x'});
    
    map = TensorMap();
    map.fromTbl = Gtbl;
    map.toTbl = Ltbl;
    map.mergefds = {'x', 'y', 'z'};

    cells = map.getDispatchInd();
    
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
