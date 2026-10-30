classdef PlatiniumKineticFunctions

    properties

        k1p
        k1m
        k2p
        k2m
        k3p
        k3m

        chargeTransferCoefficient

        th % threshold constant used for the regularization of the square root.
           % not used here
        
    end
    
    methods

        function kf = PlatiniumKineticFunctions(parameters)
            
            kf = setupParameters(kf, parameters);

        end

        function rate = computeRate(kf, eta, aH2O, aOH, aH2, T)

            th = kf.th;
            
            R = PhysicalConstants.R;
            F = PhysicalConstants.F;

            beta_HER = 0.5;
            E_0_HER = 0;

            k1_0     = 5e-10;   % [mol/(cm^2 * s)]
            k1_neg_0 = 6e-9;    % [mol/(cm^2 * s)]
            k2_0     = 1.5e-11; % [mol/(cm^2 * s)]
            k3_0     = 3e-10;   % [mol/(cm^2 * s)]

            change_unit = true;
            if change_unit
                k1_0     = 1e4*k1_0     ; 
                k1_neg_0 = 1e4*k1_neg_0 ; 
                k2_0     = 1e4*k2_0     ; 
                k3_0     = 1e4*k3_0     ; 
            end

            k2_neg_0 = k1_0*k2_0/k1_neg_0; % paper: 4.2e-12 [mol/(cm^2 * s)]
            k3_neg_0 = (k1_0^2)*k3_0/(k1_neg_0^2); % paper: 6.2e-12 [mol/(cm^2 * s)]

            k1     = k1_0.*exp(-beta_HER.*F.*(eta)./(R.*T)).*aH2O;
            k1_neg = k1_neg_0.*exp((1 - beta_HER).*F.*(eta)./(R.*T)).*aOH;
            k2     = k2_0.*exp(-beta_HER.*F.*(eta)./(R.*T)).*aH2O;
            k2_neg = k2_neg_0.*exp((1 - beta_HER).*F.*(eta)./(R.*T)).*aOH.*aH2;
            k3     = k3_0;
            k3_neg = k3_neg_0.*aH2;

            a_theta = 2*(k3 - k3_neg);
            b_theta = k1 + k1_neg + k2 + k2_neg +4*k3_neg;
            c_theta = -(k1 + k2_neg + 2*k3_neg);

            % Compute the square root, we do not use threshold
            x = b_theta.^2 - 4*a_theta.*c_theta;
            sqr = x; % quick way to create AD with same size
            ind = (x <= 0);
            sqr(ind)  = 0*x(ind);
            sqr(~ind) = (x(~ind)).^0.5;

            % Original expression is
            % theta_HER_full = (-b_theta + y)./(2*a_theta);
            % which is equivalent to
            theta_HER_full = -2*c_theta./(sqr + b_theta);

            debug = false;

            if debug

                fprintf('\n');
                kf.printstat(a_theta, 'a');
                kf.printstat(b_theta, 'b');
                kf.printstat(c_theta, 'c');
                kf.printstat(theta_HER_full, 'theta');
                
            end
            
            theta_HER_full = min(1, theta_HER_full);
            theta_HER_full = max(0, theta_HER_full);
            
            v1_full = k1.*(1 - theta_HER_full) - k1_neg.*theta_HER_full;
            v2_full = k2.*theta_HER_full - k2_neg.*(1 - theta_HER_full);

            rate = -F*(v1_full + v2_full);

        end
        
        function rate = computeIonomerRate(kf, eta, aH2O, aOH, aH2, T)

            rate = kf.computeRate(eta, aH2O, aOH, aH2, T);

        end
        
        function rate = computeElyteRate(kf, eta, aH2O, aOH, aH2, T)
            
            rate = kf.computeRate(eta, aH2O, aOH, aH2, T);
            
        end
        
    end

    methods (Static)
        
        function printstat(x, name)
            x = value(x);
            fprintf('%s, min : %g, max: %g, mean: %g\n', name, min(x), max(x), mean(x));
        end
        
    end
    
end
