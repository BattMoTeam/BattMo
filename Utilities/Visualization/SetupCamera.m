classdef SetupCamera

    properties

        G
        z
        cameraTarget
        cameraDistance = 0.1 % distance from target
        azimuthalAngle = 80 % azimuthal angle angle (phi)
        polarAngle     = 45 % polar angle angle theta (theta)
        viewAngle      = 60 % camera view angle
    end

    methods

        function cam = SetupCamera(G)
            cam.G;
            cam.z = max(G.faces.centroids(:, 3));
            cam.cameraTarget = 0.5*[0; 0; cam.z];
        end

        function do(cam)

            target    = cam.cameraTarget;
            camdist   = cam.cameraDistance;
            theta     = cam.polarAngle*(2*pi)/360;
            phi       = cam.azimuthalAngle*(2*pi)/360;
            viewangle = cam.viewAngle;

            Mtheta = [ [ cos(theta) , sin(theta), 0];
                       [ -sin(theta), cos(theta), 0];
                       [ 0, 0, 1]];
            Mphi = [[cos(phi)  , 0, sin(phi) ];
                    [0         , 1, 0        ];
                    [-sin(phi), 0, cos(phi) ]];
            M = Mtheta*Mphi;
            d = M*[0; 0; 1];
            pos = target + camdist*d;
            if phi == 0
                uv = [0; theta; 0];
            else
                uv = [0; 0; 1];
                uv = uv - (uv'*d)/(d'*d)*d;
            end

            set(gca, 'ZDir', 'normal')
            set(gca, 'CameraPosition', pos   , ...
                     'CameraTarget'  , target, ...
                     'CameraUpVector', uv)
            set(gca, 'CameraViewAngle', viewangle);

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
