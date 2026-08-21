classdef KineticFunctions

    properties

        k1p
        k1m
        k2p
        k2m
        k3p
        k3m

        chargeTransferCoefficient
        
    end
    
    methods

        function kf = KineticFunctions(parameters)
            
            kf = setupParameters(kf, parameters);
            
        end

        function rate = computeRate(kf, eta, aH2O, cOH, T)

            th = 1e-2;
            
            R = PhysicalConstants.R;
            F = PhysicalConstants.F;

            beta = kf.chargeTransferCoefficient;
            
            fp = exp(-beta.*F.*eta./(R.*T));
            fm = exp(-(1 - beta).*F.*eta./(R.*T));

            k1p_actif = kf.k1p.*aH2O.*fp;
            k1m_actif = kf.k1m.*(cOH/1000).*fm;
            k2p_actif = kf.k2p.*aH2O.*fp;
            k2m_actif = kf.k2m.*(cOH/1000).*fm;
            k3m_actif = kf.k3m; % we do not include the H2 activity
            k3p_actif = kf.k3p; % we do not include the H2 activity
            
            x = (k1p_actif + k1m_actif + k2p_actif);

            lambda = (-x + regularizedSqrt(x.^2 + 8*k1p_actif.*k3m_actif, th))./(4*k3p_actif);

            rate = -F*(k1p_actif.*(1 - lambda) - k1m_actif.*lambda + ...
                      k2p_actif.*lambda - k2m_actif.*(1 - lambda));

            r = rate;
            if isa(r, 'ADI')
                r = combineEquations(r);
                if any(isnan(r.val)) || nnz(any(isnan(r.jac{1}))) > 0
                    % keyboard
                end
            end
        end
        
        function rate = computeIonomerRate(kf, eta, aH2O, cOH, T)

            rate = kf.computeRate(eta, aH2O, cOH, T);

        end
        
        function rate = computeElyteRate(kf, eta, aH2O, cOH, T)
            
            rate = kf.computeRate(eta, aH2O, cOH, T);
            
        end
        
    end
    
end
