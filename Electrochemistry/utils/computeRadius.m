function R_lith = computeRadius(cAverage, cmax, R_delith)

        
            %theta100  = SwellingMaterialmodel.Interface.theta100;
            %theta0    = SwellingMaterialmodel.Interface.theta0;

            molarVolumeSi = 1.2e-05;
            molarVolumeLi = 9e-06;

            % We cannot anymore express <x> as a ratio of between the current and maximal concentrations of Li because the radius of the 
            % particle is changing. We have to express it as a ratio of matter quantities N/Nmax. The expression above is from a
            % rearrangement of equations 11 and 14 (can seem complex but can be done manually easily) in 'Analysis of Lithium Insertion/Deinsertion in a Silicon Electrode
            % Particle at Room Temperature Rajeswari, Chandrasekaran,Alexandre Magasinski, Gleb Yushin, and Thomas F. Fuller'*


            % defining useful quantities
            Q = (3.75.*molarVolumeLi)./(molarVolumeSi);
            c_ratio = (cAverage./cmax);
            denom = 1 - (Q .* c_ratio)/(1+Q);

            radius = R_delith ./ (denom .^ (1/3));

            if radius < R_delith
                error('R_lith has became smaller than R_delith; impossible')
            end

            R_lith = radius;

           
            %deltaTheta = theta100 - theta0;
            %%cappingValues
            %theta = max(theta, theta0);
            %theta = min (theta, theta100);



            %R_lith = R_delith .* ((deltaTheta - Q.*theta0) ./ (deltaTheta - Q .* c_ratio)) .^ (1/3);

            %soc = (theta - theta0) ./ (theta100 - theta0);

            %R_lith = R_delith .* (1 + Q .* soc);
end