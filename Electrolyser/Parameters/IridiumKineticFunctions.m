classdef IridiumKineticFunctions

    properties

        K1
        K2
         
        effectiveChargeTransferCoefficient
        
    end
    
    methods

        function kf = IridiumKineticFunctions(parameters)
            
            kf = setupParameters(kf, parameters);
            
        end

        function rate = computeRate(kf, eta, cOH, T)

            th = 1e-2;

            
            R = PhysicalConstants.R;
            F = PhysicalConstants.F;

            alpha = kf.effectiveChargeTransferCoefficient;
            K1    = kf.K1;
            K2    = kf.K2;
            
            falpha = exp((1 - alpha).*F.*eta./(R.*T));
            f      = exp(-F.*eta./(R.*T));

            rate = (2*F*K2*((cOH/1000).^alpha).*falpha)./(1 + K1*f);

        end
        
        function rate = computeIonomerRate(kf, eta, aH2O, cOH, T)

            rate = kf.computeRate(eta, cOH, T);

        end
        
        function rate = computeElyteRate(kf, eta, aH2O, cOH, T)
            
            rate = kf.computeRate(eta, cOH, T);
            
        end
        
    end
    
end
