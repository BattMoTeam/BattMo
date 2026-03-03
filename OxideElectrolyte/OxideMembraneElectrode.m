classdef OxideMembraneElectrode < BaseModel

    properties

        constants

        T     % Temperature
        Rct   % charge transfer resistance
        Eocp  % Open circuit voltage
        muEl0 % standard value of chemical electron potential
        Keh   % Equilibrium constant for the hole-electron reaction

        pO2   % O2 pressure used to compute Eocp

        % helper structure
        compinds
        
    end

    methods

        function model = OxideMembraneElectrode(inputparams)

            model = model@BaseModel();

            fdnames = {'T'     , ...
                       'pO2'   , ...
                       'muEl0', ...
                       'Keh'   , ...
                       'Rct'};
            model = dispatchParams(model, inputparams, fdnames);

            model.constants = PhysicalConstants();

            c = model.constants;

            Eocp = c.R*model.T/c.F*log(model.pO2);

            model.Eocp = Eocp;

            compinds.ch = 1;
            compinds.ce = 2;

            model.compinds = compinds;

        end

        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@BaseModel(model);

            varnames = {};

            % Electrode electrostatic potential
            varnames{end + 1} = 'phi';
            % Electrode electromotive potential
            varnames{end + 1} = 'pi';
            % Electrode electronic chemical potential
            varnames{end + 1} = 'E';
            % Electron concentration
            varnames{end + 1} = 'ce';
            % Hole concentration
            varnames{end + 1} = 'ch';
            % Logarithmic concentrations
            varnames{end + 1} = VarName({}, 'logcs', 2);
            % Equation for pi - phi
            varnames{end + 1} = 'currentVoltageEquation';
            % Equilibrium equation for hole-electron reaction
            varnames{end + 1} = 'equilibriumEquation';
            % Electronic flux at electrode (positif is flux leaving electrode)
            varnames{end + 1} = 'jEl';
            % ionix flux at electrode (positif is flux leaving electrode)
            varnames{end + 1} = 'jO2m';
            % Charge flux at electrode (positif is flux leaving electrode)
            varnames{end + 1} = 'j';
            % Charge conservation equation (that is j - jEl - jHp = 0)
            varnames{end + 1} = 'chargeCons';
            % Definition equation for jEl
            varnames{end + 1} = 'jElEquation';
            % Definition equation for jO2m
            varnames{end + 1} = 'jO2mEquation';
            % Electronic concentration definition equation
            varnames{end + 1} = 'electronConcentrationEquation';
            model = model.registerVarNames(varnames);

            fn = @OxideMembraneElectrode.updateConcentrations;
            inputnames = {VarName({}, 'logcs', 2)};
            model = model.registerPropFunction({'ce', fn, inputnames});
            model = model.registerPropFunction({'ch', fn, inputnames});

            fn = @OxideMembraneElectrode.updateChargeCons;
            inputnames = {'j', 'jEl', 'jO2m'};
            model = model.registerPropFunction({'chargeCons', fn, inputnames});

            fn = @OxideMembraneElectrode.updateE;
            inputnames = {'phi', 'pi'};
            model = model.registerPropFunction({'E', fn, inputnames});

            fn = @OxideMembraneElectrode.updateCurrentVoltageEquation;
            inputnames = {'j', 'E'};
            model = model.registerPropFunction({'currentVoltageEquation', fn, inputnames});

            fn = @OxideMembraneElectrode.updateElConcEquation;
            inputnames = {'E', 'ce'};
            model = model.registerPropFunction({'electronConcentrationEquation', fn, inputnames});

            fn = @OxideMembraneElectrode.updateEquilibriumReaction;
            inputnames = {'ch', 'ce'};
            model = model.registerPropFunction({'equilibriumEquation', fn, inputnames});

        end

        function state = updateConcentrations(model, state)

            compinds = model.compinds;

            logcs = state.logcs;

            state.ch = exp(logcs{compinds.ch});
            state.ce = exp(logcs{compinds.ce});

        end


        function state = updateChargeCons(model, state)

            j   = state.j;
            jEl = state.jEl;
            jO2m = state.jO2m;

            state.chargeCons = j - (jEl - 2*jO2m);

        end

        function state = updateE(model, state)

            pi   = state.pi;
            phi  = state.phi;

            state.E = pi - phi;

        end

        function state = updateCurrentVoltageEquation(model, state)

            Eocp = model.Eocp;
            Rct  = model.Rct;

            E = state.E;
            j = state.j;

            cveq = (E - Eocp)/Rct - j;

            state.currentVoltageEquation = cveq;

        end

        function state = updateElConcEquation(model, state)

            mu0 = model.muEl0;
            c   = model.constants;
            T   = model.T;

            ce = state.ce;
            E  = state.E;

            eceq = E + (1/c.F)*(mu0 + (c.R*T)*log(ce));

            state.electronConcentrationEquation = eceq;

        end


        function state = updateEquilibriumReaction(model, state)

            Keh = model.Keh;

            ch = state.ch;
            ce = state.ce;

            state.equilibriumEquation = ch*ce/Keh - 1;

        end
    end

end
