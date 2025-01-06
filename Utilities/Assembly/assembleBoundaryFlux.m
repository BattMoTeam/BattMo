function [boundaryFlux, extra] = assembleBoundaryFlux(model, potential, boundary_potential, fluxCoefficient, boundaryFaces)
% Compute and returns the flux a boundary. 
%
% The boundary is described by coupterm which is an instance of couplingTerm
%
% The sign of the flux is oriented in the exterior direction
%
% The flux corresponds to the discretization of  $ - \int_{face} D\nabla(p)\cdot n ds$, where
% - D                  : flux coefficient
% - p                  : potential
% - \nable             : gradient operator
% - n                  : normal pointing towards the exterior
% - \int_{face} ... ds : integrale over a boundary face
% Note the minus sign (same convention as in assembleFlux)
%
% Additional output in variable extra is returned, to be used in assembleBoundarySource

    bcval = boundary_potential;
    
    % Retrieve the half-transmissibilities at the boundary faces
    [halftrans, boundaryCells, sgn] = model.G.getBcTrans(boundaryFaces);

    val = potential(boundaryCells);
    
    % Compute boundary flux
    boundaryFlux = fluxCoefficient.*halftrans.*(bcval - val);

    if nargout > 1
        extra = struct('sgn', sgn, ...
                       'boundaryCells', boundaryCells);
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
