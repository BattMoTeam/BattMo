classdef SetupCamera

    properties

        G
        cameraTarget
        cameraDistance
        azimuthalAngle = 80 % azimuthal angle angle
        polarAngle     = 45 % polar angle angle
        viewAngle      = 60 % camera view angle
    end

    methods

        function cam = SetupCamera(G)
            cam.G;
            cam.cameraTarget = [0; 0; 0];
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
                    [-sin(phi)], 0, cos(phi) ];
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

