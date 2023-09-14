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
            inputnames = {'ch', 'ce', 'phi'};
            model = model.registerPropFunction({'jEl', fn, inputnames});

            fn = @OxideMembraneElectrolyte.updateMassConsO2;
            inputnames = {'sourceO2', 'jO2'};
            model = model.registerPropFunction({'massConsO2', fn, inputnames});

            fn = @OxideMembraneElectrolyte.updateChargeConsEl;
            inputnames = {'sourceEl', 'jEl'};
            model = model.registerPropFunction({'chargeConsEl', fn, inputnames});

        end


        function state = updateO2Flux(model, state)

            sigmaO2 = model.sigmaO2;
            op      = model.operators;

            phi = state.phi;

            state.jO2 = assembleHomogeneousFlux(model, phi, sigmaO2);

        end


        function state = updateElFlux(model, state)

            c  = model.constants;
            Dh = model.Dh;
            De = model.De;
            
            phi = state.phi;
            ce  = state.ce;
            ch  = state.ch;
            
            jch = assembleHomogenousFlux(model, ch, c.F*Dh);
            jce = assembleHomogenousFlux(model, ce, c.F*De);

            effSigma = (c.F)^2/(c.R*c.T)*(Dh*ch + De*ce);

            jphi = assembleFlux(model, phi, effSigma);

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
