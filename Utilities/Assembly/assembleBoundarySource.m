function [jExternal, jFaceExternal] = assembleBoundarySource(model, potential, boundary_potential, fluxCoefficient, coupterm)
% Returns a the source term (in fact a "sink") which correponds to the flux that leaves the domain at the boundary
%
% See function assembleBoundaryFlux
%

    [boundaryFlux, extra] = assembleBoundaryFlux(model, potential, boundary_potential, fluxCoefficient, coupterm);
    
    faces = coupterm.couplingfaces;

    G = model.G;
    nc = G.topology.cells.num;
    zeroFaceAD = model.AutoDiffBackend.convertToAD(zeros(nc, 1), potential);    
    jExternal = zeroFaceAD;
    jExternal = subsetPlus(jExternal, boundaryFlux, extra.cells);

    if nargout > 1
        nf = G.topology.faces.num;
        zeroFaceAD = model.AutoDiffBackend.convertToAD(zeros(nf, 1), potential);
        jFaceExternal = zeroFaceAD;
        jFaceExternal = subsasgnAD(jFaceExternal, faces, - extra.sgn.*boundaryFlux);
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
