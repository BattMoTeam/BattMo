classdef IonomerMembrane < ElectronicComponent
    
    properties
        
        liquidVolumeFraction
        
        H2O % with fields
        % H2O.c0 :  initial reference concentration
        % H2O.D :  diffusion coefficient for water
        
        OH % with fields
        % OH.xi : OH occupation
        % OH.z :  
        % OH.t :
    
    end
    
    methods
        
        function model = IonomerMembrane(paramobj)
            
            model = model@ElectronicComponent(paramobj);
        end
        
        
        function state = updateH2Oc(model, state)
            H2Oceps = state.H2Oceps;
            
            eps = model.liquidVolumeFraction;
            
            state.H2Oc = H2Oceps./eps;
        end

        
        function state = updateH2Oactivity(model, state)
            H2Oc = state.H2Oc;

            H2Oc0 = model.H2O.c0;
            
            state.H2Oactivity = H2Oc./H2Oc0;
        end
        
        
        function state = updateKappaEff(model, state)
            T = state.T;
            aw = state.H2Oactivity;
            
            eps = model.liquidVolumeFraction;
            
            aw(aw > 1) = 1;
            
            kappa = (0.1334 - 3.882e-4.*T + (0.01148.*T - 3.909).*aw ...
                     - (0.06690.*T - 23.01).*aw.^2 + (0.1227.*T - 42.61).*aw.^3 ...
                     - (0.06021.*T - 21.80).*aw.^4) .* 20;
            
            state.kappaeff = kappa.*eps.^1;
        end
        
        
        function state = updatedmudc(model, state)
            T = state.T;
            c = state.H2Oc;
            
            R = model.constants.R;
            
            state.H2Odmudc = R.*T./c;
        end
        
        
        function state = updatejchem(model, state)
            kappaeff = state.kappaeff;
            dmudc    = state.H2Odmudc;
            c        = state.H2Oc;
            
            xi = model.OH.xi;
            z  = model.OH.z;
            F  = model.constants.F;
            
            fluxcoef = dmudc.*kappaeff.*xi./(z.*F);
            
            state.jchem = assembleFlux(model, c, fluxcoef);
        end
        
        
        function state = updatej(model, state)
            kappaeff = state.kappaeff;
            phi      = state.phi;
            jchem    = state.jchem;
            
            fluxcoef = dmudc.*kappaeff.*xi./(z.*F);
     
            j = assembleFlux(model, phi, kappeff);
            
            state.j = j + jchem;
        end

        function state = updateH2OdiffFlux(model, state)
            c = state.H2Oc;
            
            D = model.H2O.D;
            eps = model.liquidVolumeFraction;
            
            Deff = D.*eps.^1;
            
            state.H2OdiffFlux = assembleFlux(model, c, Deff);
        end
    
    
        function state = updateH2OmigFlux(model, state)
            j = state.j;
            
            xi = model.OH.xi;
            z  = model.OH.z;
            t  = model.OH.t;
            F  = model.constants.F;
            
            state.H2OmigFlux = xi.*t/(z.*F).*j;
        end
    
        
        function state = assembleH2OmassCon(model, state)
            
            diffFlux = state.H2OdiffFlux;
            migFlux  = state.H2OmigFlux;
            accum    = state.H2Oaccum;
            source   = state.H2OSource;
            
            flux = difflux + migFlux;
            bcsource  = 0;
            
            state.H2OmassCons = assembleConservationEquation(model, flux, bcsource, source, accum)
        end
        
    
    end
    
end

