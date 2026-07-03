classdef OxideMembraneElectrolyte < BaseModel
%
% ref1 : Jacobsen and Mogensen (2008) : The course of oxygen partial pressure and electric potentials across oxide electrolyte cell
%
    properties

        % Temperature
        T
        % Structure with physical constants
        constants
        % Diffusion constant for hole
        Dh
        % Diffusion constant for electron
        De
        % O2- conductivity
        sigmaO2
        % Equilibrium Constant
        Keh
        % standard value of chemical electron potential
        muEl0
        
        % Helper
        compinds

    end

    methods

        function model = OxideMembraneElectrolyte(inputparams)

            model = model@BaseModel();

            fdnames = {'G'      , ...
                       'T'      , ...
                       'Dh'     , ...
                       'De'     , ...
                       'sigmaO2', ...
                       'Keh'    , ...
                       'muEl0'};

            model = dispatchParams(model, inputparams, fdnames);

            model.operators = localSetupOperators(model.G);

            model.constants = PhysicalConstants();

            con = model.constants;

            compinds.ch = 1;
            compinds.ce = 2;

            model.compinds = compinds;

        end

        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@BaseModel(model);

            varnames = {};

            % Electrostatic potential phi
            varnames{end + 1} = 'phi';
            % Concentration electron
            varnames{end + 1} = 'ce';
            % Concentration hole
            varnames{end + 1} = 'ch';
            % Logarithmic concentrations
            varnames{end + 1} = VarName({}, 'logcs', 2);
            % Electronic Flux in Ampere (multiplied with faceArea)
            varnames{end + 1} = 'jEl';
            % Ionic Flux in Ampere (multiplied with faceArea)
            varnames{end + 1} = 'jO2m';
            % Coefficient in front of grad(phi) in assembly of jEl
            varnames{end + 1} = 'gradPhiCoef';
            % Coefficient in front of grad(ce) and grad(ch) assembly of jEl
            varnames{end + 1} = VarName({}, 'gradConcCoefs', 2);
            % H+ source term
            varnames{end + 1} = 'sourceO2';
            % regularization parameter
            varnames{end + 1} = 'alpha';
            % Electronic source term
            varnames{end + 1} = 'sourceEl';
            % O2 mass conservation (we measure the mass in Coulomb, hence "massConsHp")
            varnames{end + 1} = 'massConsO2';
            % Charge conservation
            varnames{end + 1} = 'chargeConsEl';
            % Equilibrium equation for hole-electron reaction
            varnames{end + 1} = 'equilibriumEquation';
            % electromotive potential (not used in assembly of equations)
            varnames{end + 1} = 'pi';
            
            model = model.registerVarNames(varnames);

            % The variable pi is not used for the assembly of the residual equations. We register it as such.
            varname = 'pi';
            model = model.registerExtraVarName(varname);
            
            fn = @OxideMembraneElectrode.updateConcentrations;
            inputnames = {VarName({}, 'logcs', 2)};
            model = model.registerPropFunction({'ce', fn, inputnames});
            model = model.registerPropFunction({'ch', fn, inputnames});

            fn = @OxideMembraneElectrolyte.updateO2Flux;
            inputnames = {'phi'};
            model = model.registerPropFunction({'jO2m', fn, inputnames});

            fn = @OxideMembraneElectrolyte.updateEquilibriumEquation;
            inputnames = {'ch', 'ce'};
            model = model.registerPropFunction({'equilibriumEquation', fn, inputnames});

            fn = @OxideMembraneElectrolyte.updateElFlux;
            inputnames = {VarName({}, 'gradConcCoefs', 2), ...
                          'gradPhiCoef'                  , ...
                          'phi'                          , ...
                          'ce'                           , ...
                          'ch'};
            model = model.registerPropFunction({'jEl', fn, inputnames});

            fn = @OxideMembraneElectrolyte.updateMassConsO2;
            inputnames = {'sourceO2', 'jO2m'};
            model = model.registerPropFunction({'massConsO2', fn, inputnames});

            fn = @OxideMembraneElectrolyte.updateChargeConsEl;
            inputnames = {'sourceEl', 'jEl'};
            model = model.registerPropFunction({'chargeConsEl', fn, inputnames});

            fn = @OxideMembraneElectrolyte.updateGradCoefs;
            inputnames = {'ce', 'ch', 'alpha'};
            model = model.registerPropFunction({'gradPhiCoef', fn, inputnames});
            model = model.registerPropFunction({VarName({}, 'gradConcCoefs', 2), fn, inputnames});

            fn = @OxideMembraneElectrolyte.updatePi;
            inputnames = {'ce', 'phi'};
            model = model.registerPropFunction({'pi', fn, inputnames});

        end

        function state = updateConcentrations(model, state)

            compinds = model.compinds;

            logcs = state.logcs;

            state.ch = exp(logcs{compinds.ch});
            state.ce = exp(logcs{compinds.ce});

        end

        function state = updateGradCoefs(model, state)

            Dh    = model.Dh;
            De    = model.De;
            T     = model.T;
            c     = model.constants;
            cinds = model.compinds;

            ce    = state.ce;
            ch    = state.ch;
            alpha = state.alpha;

            chr = exp(alpha*log(ch));
            cer = exp(alpha*log(ce));

            gConcCs{cinds.ch} = c.F*Dh*(chr./ch);
            gConcCs{cinds.ce} = c.F*De*(cer./ce);

            gPhiC = (c.F)^2/(c.R*T)*(Dh*chr + De*cer);

            state.gradConcCoefs = gConcCs;
            state.gradPhiCoef   = gPhiC;

        end


        function state = updateO2Flux(model, state)

            sigmaO2 = model.sigmaO2;
            op      = model.operators;

            phi = state.phi;

            % Note minus sign in front
            state.jO2m = - assembleHomogeneousFlux(model, phi, sigmaO2);

        end


        function state = updateElFlux(model, state)
            % Equation (17) in ref1
            cinds = model.compinds;
            
            phi     = state.phi;
            ce      = state.ce;
            ch      = state.ch;
            gPhiC   = state.gradPhiCoef;
            gConcCs = state.gradConcCoefs;

            jch  = assembleFlux(model, ch, gConcCs{cinds.ch});
            jce  = assembleFlux(model, ce, gConcCs{cinds.ce});
            jphi = assembleFlux(model, phi, gPhiC);

            state.jEl = jch - jce + jphi;

        end


        function state = updateChargeConsEl(model, state)

            op = model.operators;

            sourceEl = state.sourceEl;
            jEl      = state.jEl;

            state.chargeConsEl =  op.Div(jEl) - sourceEl;

        end

        function state = updateMassConsO2(model, state)

            op = model.operators;

            sourceO2 = state.sourceO2;
            jO2m     = state.jO2m;

            state.massConsO2 =  op.Div(jO2m) - sourceO2;

        end

        function state = updateEquilibriumEquation(model, state)

            Keh = model.Keh;

            ch = state.ch;
            ce = state.ce;

            state.equilibriumEquation = ch.*ce/Keh - 1;

        end

        function state = updatePi(model, state)

            c   = model.constants;
            T   = model.T;
            mu0 = model.muEl0;
            
            phi = state.phi;
            ce  = state.ce;
            
            state.pi = 1/c.F*(- mu0 - c.R*T*log(ce)) + phi;
            
        end
        
    end

end
