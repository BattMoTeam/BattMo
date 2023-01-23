classdef IonomerMembrane < ElectronicComponent
    
    properties
        
        volumeFraction
        
        H2O % with fields
        % H2O.c0 : Reference concentration
        % H2O.D  : diffusion coefficient for water
        % H2O.V0 : Molar mass (needed for function groupHydration which is only needed in setup of initial condition and not for assembly)
        
        OH % with fields
        % OH.xi : OH occupation
        % OH.z  : charge number 
        % OH.t  : transference number

        cT % Total concentration of charged groups

        V % molar volume (needed for function groupHydration which is only needed in setup of initial condition and not for assembly)
        
    end
    
    methods
        
        function model = IonomerMembrane(paramobj)

            paramobj.use_thermal = false;
            model = model@ElectronicComponent(paramobj);

            fdnames = {'volumeFraction', ...
                       'H2O'           , ...
                       'OH'            , ...
                       'cT'            , ...
                       'V'};
            model = dispatchParams(model, paramobj, fdnames);

            model.constants = PhysicalConstants();
            
        end

        
        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@ElectronicComponent(model);

            varnames = {};
            
            % Water concentration (per total volume)
            varnames{end + 1} = 'H2Oceps';
            % Water concentration
            varnames{end + 1} = 'H2Oc';
            % OH concentration
            varnames{end + 1} = 'cOH';
            
            %% Component properties
            
            % H2O activity
            varnames{end + 1} = 'H2Oa'; 
            % Chemical potential derivative
            varnames{end + 1} = 'H2Odmudc';
            
            %% Flux variables
            
            % Chemical flux
            varnames{end + 1} = 'jchem';
            % H2O diffusion flux
            varnames{end + 1} = 'H2OdiffFlux';
            % H2O migration flux
            varnames{end + 1} = 'H2OmigFlux';
            
            %% Source terms
            
            % H2O source [mol/s] (one value per cell)
            varnames{end + 1} = 'H2OSource';
            % OH source [mol/s] (one value per cell)
            varnames{end + 1} = 'OHSource';
        
            %% Accumulation term
            
            % H2O accumulation term
            varnames{end + 1} = 'H2Oaccum';
            
            %% Mass conservation equation
            
            % H2O mass conservation equation
            varnames{end + 1} = 'H2OmassCons';

            model = model.registerVarNames(varnames);

            % update water concentration
            fn = @() IonomerMembrane.updateH2Oc;
            inputnames = {'H2Oceps'};
            model = model.registerPropFunction({'H2Oc', fn, inputnames});
            
            % update water activity
            fn = @() IonomerMembrane.updateH2Oactivity;
            inputnames = {'H2Oc'};
            model = model.registerPropFunction({'H2Oa', fn, inputnames});
            
            % update effective conductivity
            fn = @() IonomerMembrane.updateConductivity;
            inputnames = {'T', 'H2Oa'};
            model = model.registerPropFunction({'conductivity', fn, inputnames});

            % update chemical potential derivative
            fn = @() IonomerMembrane.updatedmudc;
            inputnames = {'T', 'H2Oc'};
            model = model.registerPropFunction({'H2Odmudc', fn, inputnames});
            
            % update chemical flux
            fn = @() IonomerMembrane.updatejchem;
            inputnames = {'conductivity', 'H2Odmudc', 'H2Oc'};
            model = model.registerPropFunction({'jchem', fn, inputnames});            

            % update electronic flux
            fn = @() IonomerMembrane.updateCurrent;
            inputnames = {'conductivity', 'phi', 'jchem'};
            model = model.registerPropFunction({'j', fn, inputnames});            

            % update H2O diffusion flux
            fn = @() IonomerMembrane.updateH2OdiffFlux;
            inputnames = {'H2Oc'};
            model = model.registerPropFunction({'H2OdiffFlux', fn, inputnames});
            
            % update H2O migration flux
            fn = @() IonomerMembrane.updateH2OmigFlux;
            inputnames = {'j'};
            model = model.registerPropFunction({'H2OmigFlux', fn, inputnames});
            
            % update electronic source
            fn = @() IonomerMembrane.updateEsource;
            inputnames = {'OHSource'};
            model = model.registerPropFunction({'eSource', fn, inputnames});

            % update electronic source
            fn = @() IonomerMembrane.assembleH2Oaccum;
            functionCallSetupFn = @(propfunction) PropFunction.accumFuncCallSetupFn(propfunction);
            fn = {fn, functionCallSetupFn};
            inputnames = {'H2Oceps'};
            model = model.registerPropFunction({'H2Oaccum', fn, inputnames});
                        
            % assemble H2O mass conservation equation
            fn = @() IonomerMembrane.assembleH2OmassCons;
            inputnames = {'H2OdiffFlux', 'H2OmigFlux', 'H2Oaccum', 'H2OSource'};
            model = model.registerPropFunction({'H2OmassCons', fn, inputnames});

            % the OH concentration is constant
            fn = @() IonomerMembrane.setupOHconcentration;
            inputnames = {};
            model = model.registerPropFunction({'cOH', fn, inputnames});
            
        end

        function state = updateH2Oc(model, state)

            H2Oceps = state.H2Oceps;
            
            vf = model.volumeFraction;
            
            state.H2Oc = H2Oceps./vf;
            
        end

        
        function state = setupOHconcentration(model, state)

            state.cOH = model.cT;
            
        end
        
        
        function state = updateH2Oactivity(model, state)

            H2Oc = state.H2Oc;
            H2Oc0 = model.H2O.c0;
            
            state.H2Oa = H2Oc./H2Oc0;
            
        end
        
        
        function state = updateConductivity(model, state)
            
            T = state.T;
            a = state.H2Oa;
            
            vf = model.volumeFraction;
            
            a(a > 1) = 1;
            
            kappa = (0.1334                       ...
                     - 3.882e-4.*T                ...
                     + (0.01148.*T - 3.909).*a    ...
                     - (0.06690.*T - 23.01).*a.^2 ...
                     + (0.1227.*T - 42.61).*a.^3  ...
                     - (0.06021.*T - 21.80).*a.^4) .* 20;
            
            state.conductivity = kappa.*vf.^1;
            
        end
        
        
        function state = updatedmudc(model, state)
            
            T = state.T;
            c = state.H2Oc;
            
            R = model.constants.R;
            
            state.H2Odmudc = R.*T./c;
            
        end
        
        
        function state = updatejchem(model, state)
            
            conductivity = state.conductivity;
            dmudc        = state.H2Odmudc;
            c            = state.H2Oc;
            
            xi = model.OH.xi;
            z  = model.OH.z;
            F  = model.constants.F;
            
            %% TODO : check expression (z^2F^2 in paper)
            fluxcoef = dmudc.*conductivity.*xi./(z.*F);
            
            state.jchem = assembleFlux(model, c, fluxcoef);
        end
        
        
        function state = updateCurrent(model, state)
            
            conductivity = state.conductivity;
            phi          = state.phi;
            jchem        = state.jchem;
            
            jelec = assembleFlux(model, phi, conductivity);
            
            state.j = jelec + jchem;
            
        end

        function state = updateH2OdiffFlux(model, state)

            D  = model.H2O.D;
            vf = model.volumeFraction;

            c = state.H2Oc;
            
            Deff = D.*vf.^1;
            
            state.H2OdiffFlux = assembleFlux(model, c, Deff);
            
        end
    
    
        function state = updateH2OmigFlux(model, state)

            j = state.j;
            
            xi = model.OH.xi;
            z  = model.OH.z;
            t  = model.OH.t;
            F  = model.constants.F;
            
            %% todo : check expression, not same as in paper
            state.H2OmigFlux = xi.*t/(z.*F).*j;
            
        end
    
        function state = assembleH2OmassCons(model, state)
            
            diffFlux = state.H2OdiffFlux;
            migFlux  = state.H2OmigFlux;
            accum    = state.H2Oaccum;
            source   = state.H2OSource;
            
            flux = difflux + migFlux;
            bcsource  = 0;
            
            state.H2OmassCons = assembleConservationEquation(model, flux, bcsource, source, accum)
        end
        
        function state = updateEsource(model, state)

            F = model.constants.F;
            z = model.OH.z;
            
            state.eSource = z*F*state.OHSource;
            
        end

        function state = assembleH2Oaccum(model, state, state0, dt)

            c  = state.H2Oceps;
            c0 = state0.H2Oceps;
            
            vols = model.G.cells.volumes;
            
            state.H2Oaccum = 1/dt*vols.*(c - c0);
            
        end


        
    end

    methods(Static)

        function cH2O = groupHydration(model, aw, T)

            V = model.V;
            V0 = model.H2O.V0;
            
            aw(aw > 1) = 1;
            %   lambda corresponds to the number of water molecules per charged
            %   group, in the membrane and ionomer phase, given the membrane water
            %   activity and the temperature. The default polynomial relationship is based on
            %   Jiao et al., Int J Hydrogen Energy 2014, doi:10.1016/j.ijhydene.2014.01.180
            %   using only the water activity < 1 portion of their equation.
            lambda = ((-0.6.*aw.^3 + 0.85.*aw.^2 - 0.2.*aw + 0.153) .* (T - 313) ...
                      + 39.*aw.^3 - 47.7.*aw.^2 + 23.4.*aw + 0.117).*29./18.877;

            cH2O = lambda./(V0.*lambda + V);
            
        end
        
    end
    
end

