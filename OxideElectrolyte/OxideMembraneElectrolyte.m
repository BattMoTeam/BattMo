classdef OxideMembraneElectrolyte < BaseModel

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

        % Helper
        compinds

    end

    methods

        function model = OxideMembraneElectrolyte(paramobj)

            model = model@BaseModel();

            fdnames = {'G'      , ...
                       'T'      , ...
                       'Dh'     , ...
                       'De'     , ...
                       'sigmaO2', ...
                       'Keh'};

            model = dispatchParams(model, paramobj, fdnames);

            model.operators = localSetupOperators(model.G);

            model.constants = PhysicalConstants();

            con = model.constants;

            compinds.ch = 1;
            compinds.ce = 2;

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
            % Electronic Flux
            varnames{end + 1} = 'jEl';
            % Ionic Flux
            varnames{end + 1} = 'jO2';
            % Coefficient in front of grad(phi) in assembly of jEl
            varnames{end + 1} = 'gradPhiCoef';
            % Coefficient in front of grad(ce) and grad(ch) assembly of jEl
            varnames{end + 1} = VarName({}, 'gradConcCoefs', 2);
            % H+ source term
            varnames{end + 1} = 'sourceO2';
            % Electronic source term
            varnames{end + 1} = 'sourceEl';
            % O2 mass conservation (we measure the mass in Coulomb, hence "massConsHp")
            varnames{end + 1} = 'massConsO2';
            % Charge conservation
            varnames{end + 1} = 'chargeConsEl';
            % Equilibrium equation for hole-electron reaction
            varnames{end + 1} = 'equilibriumEquation';

            model = model.registerVarNames(varnames);

            fn = @OxideMembraneElectrolyte.updateO2Flux;
            inputnames = {'phi'};
            model = model.registerPropFunction({'jO2', fn, inputnames});

            fn = @OxideMembraneElectrolyte.updateEquilEquation;
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
            inputnames = {'sourceO2', 'jO2'};
            model = model.registerPropFunction({'massConsO2', fn, inputnames});

            fn = @OxideMembraneElectrolyte.updateChargeConsEl;
            inputnames = {'sourceEl', 'jEl'};
            model = model.registerPropFunction({'chargeConsEl', fn, inputnames});

            fn = @OxideMembraneElectrolyte.updateGradCoefs;
            inputnames = {'ce', 'ch'};
            model = model.registerPropFunction({'gradPhiCoef', fn, inputnames});
            model = model.registerPropFunction({VarName({}, 'gradConcCoefs', 2), fn, inputnames});

        end

        function state = updateGradCoefs(model, state)

            Dh    = model.Dh;
            De    = model.De;
            T     = model.T;
            c     = model.constants;
            cinds = model.compinds;

            ce = state.ce;
            ch = state.ch;

            gConcCs{cinds.ch} = c.F*Dh*ch./ch;
            gConcCs{cinds.ce} = c.F*Dh*ce./ce;

            gPhiC = (c.F)^2/(c.R*c.T)*(Dh*ch + De*ce);

            state.gradConcCoefs = gConcCs;
            state.gradPhiCoef   = gPhiC;

        end


        function state = updateO2Flux(model, state)

            sigmaO2 = model.sigmaO2;
            op      = model.operators;

            phi = state.phi;

            state.jO2 = assembleHomogeneousFlux(model, phi, sigmaO2);

        end


        function state = updateElFlux(model, state)

            phi = state.phi;
            ce  = state.ce;
            ch  = state.ch;

            gPhiC   = state.(elyte).gradPhiCoef;
            gConcCs = state.(elyte).gradConcCoefs;

            jch  = assembleFlux(model, ch, gPhiCs{1});
            jce  = assembleFlux(model, ce, gPhiCs{2});
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
            jO2      = state.jO2;

            state.massConsO2 =  op.Div(jO2) - sourceO2;

        end

        function state = updateEquilibriumReaction(model, state)

            Keh = model.Keh;

            ch = state.ch;
            ce = state.ce;

            state.equilibriumEquation = ch*ce - Keh;

        end

    end

end
