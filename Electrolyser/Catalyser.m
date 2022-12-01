classdef Catalyser < BaseModel
    
    properties
        
        constants
        
        compInd   % mapping structure for component indices
        phaseInd  % mapping structure for phase indices
        liquidInd % mapping structure for component indices
        gasInd    % mapping structure for component indices

        E0 % Standard equilibrium potential,   [V]
        j0 % Exchange current density
        inmr %
        % inmr.cT   % Total concentration of charged groups
        % inmr.kxch % Rate constant for exchange between ionomer and electrolyte. [s^-1]
        % inmr.OH.z
        PF1 % should be initialized using zero
        
    end
    
    methods
        
        function model = Catalyser(paramobj)
            
        end

        function model = registerVarAndPropfuncNames(model)
            
        %% Declaration of the Dynamical Variables and Function of the model
        % (setup of varnameList and propertyFunctionList)
            
            model = registerVarAndPropfuncNames@BaseModel(model);
            
            gasElyteInd.H2O = 1;
            gasElyteInd.H2 = 2;
            ngas_elyte = 2;
            
            names = {};
            % Temperature
            varnames = {'T'};
            % Electrode potential 
            varnames{end + 1} = {'E'};
            % Electric Potential from electrolyte
            varnames{end + 1} = 'phiElyte';
            % Electric Potential from ionomer
            varnames{end + 1} = 'phiInmr';
            % concentration of proton from electrotyte
            varnames{end + 1} = 'cHElyte';
            % concentration of OH from electrotyte
            varnames{end + 1} = 'cOHElyte';
            % Electrode reactive potential 
            varnames{end + 1} = {'Er'};            
            % Water activity from electrolyte
            varnames{end + 1} = 'H2OaElyte';
            % Water activity from ionomer
            varnames{end + 1} = 'H2OaInmr';
            % Partial pressures in the electrolyte gas phase
            gasPressureElyte =  VarName({}, 'gasPressuresElyte', ngas_elyte);
            varnames{end + 1} = gasPressureElyte;

            % Potential difference on the electrolyte side
            varnames{end + 1} = {'etaElyte'};
            % Potential difference on the ionomer side
            varnames{end + 1} = {'etaInmr'};

            %% Fluxes
            
            % flux in electrolyte
            varnames{end + 1} = 'jElyte';
            % flux in ionomer
            varnames{end + 1} = 'jInmr';
            % electronic flux
            varnames{end + 1} = 'j';

            %% Rates
            
            % Ionomer ion exchange rate
            varnames{end + 1} = 'Rxch';
            % Ionomer sorption rate
            varnames{end + 1} = 'Rsorption';

            model = model.registerVarNames(varnames);
            
            % Assemble reactive potential
            fn = @() Catalyser.updateEr;
            inputnames = {'T', 'cHElyte'};
            model = model.registerPropFunction({'Er', fn, inputnames});            
            
            % Assemble reactive potential
            fn = @() Catalyser.updateEtas;
            inputnames = {'phiElyte', 'phiInmr', 'E', 'Er'};
            model = model.registerPropFunction({'eta', fn, inputnames});            
            model = model.registerPropFunction({'etaElyte', fn, inputnames});
            model = model.registerPropFunction({'etaInmr', fn, inputnames});

            % Assemble the electric fluxes
            fn = @() Catalyser.updateFluxes;
            inputnames = {'etaElyte', 'etaInmr', 'T', gasPressureElyte, 'H2OaInmr'};
            model = model.registerPropFunction({'jElyte', fn, inputnames});            
            model = model.registerPropFunction({'jInmr', fn, inputnames});            
            model = model.registerPropFunction({'j', fn, inputnames});
            
            % Assemble ionomer ion exchange rate
            fn = @() Catalyser.updateIonomerExchange;
            inputnames = {'eta', 'T', 'cOHElyte'};
            model = model.registerPropFunction({'Rxch', fn, inputnames});
            
            % Assemble sorption rate
            fn = @() Catalyser.updateSorption;
            inputnames = {'H2OaElyte', 'H2OaInmr'};
            model = model.registerPropFunction({'Rsorption', fn, inputnames});


        end
        
        function state = updateEr(model, state)

            T = state.T;
            cHElyte = state.CHElyte;
            
            E0  = model.E0;
            con = model.constants;
            F   = con.F;
            c0  = con.c0;
            R   = con.R;
            
            state.Er = E0 - R.*T./(2.*F) .* log(c0.^2 .* cHElyte.^-2);
        end
        
        function state = updateSorption(model, state)

            H2OaElyte = state.H2OaElyte;
            H2OaInmr = state.H2OaInmr;

            kML = model.inmr.kML;
            
            state.Rsorption = ionomerSorption(kML, H2OaElyte, H2OaInmr);
        end
        
        function state = updateEtas(model, state)
            
            phiElyte = state.phiElyte;
            phiInmr  = state.phiInmr;
            E        = state.E;
            Er       = state.Er;
            
            state.eta      = phiInmr - phiElyte;
            state.etaElyte = E - phiElyte - Er;
            state.etaInmr  = E - phiInmr - Er;
        end
        
        function state = updateIonomerExchange(model, state)
            
            eta      = state.eta;
            T        = state.T;
            cOHElyte = state.cOHElyte;
            
            ir = model.inmr;
            
            state.Rxch = ionomerExchange(ir.kxch, ir.cT, ir.OH.z, eta, T, cOHElyte);
        end
        
        function state = updateFluxes(model, state)
            
         etaElyte = state.etaElyte;
         etaInmr  = state.etaInmr;
         H2OaInmr = state.H2OaInmr;
         T        = state.T;
         gasp1    = state.gasPressuresElyte{1};
         gasp2    = state.gasPressuresElyte{2};
        
         ir    = model.inmr;
         Asp   = model.Asp;
         PF1   = model.PF1;
         j0    = model.j0;
         alpha = model.alpha;
         n     = model.n;
         
         X_inmr = 0.99;
         
         PF = (ir.cT/1000).*(gasp1./1e5).*(gasp2./1e5);
         jElyte = (1 - X_inmr).*Asp'.*butlerVolmer(j0, alpha, n, etaElyte, T);
         
         PF2 = H2OaInmr;
         jInmr = X_inmr.*Asp'*butlerVolmer_split(j0, PF1, PF2, alpha, n, etaInmr, T);
         
         state.jElyte = jElyte;
         state.jInmr  = jInmr;
         state.j      = jElyte + jInmr;
        end
        
        
    end

end

