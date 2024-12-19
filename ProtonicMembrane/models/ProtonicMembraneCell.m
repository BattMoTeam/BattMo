classdef ProtonicMembraneCell < BaseModel
    
    properties
        
        % Temperature
        T
        % Structure with physical constants
        constants

        Electrolyser
        GasSupply
        Interface
        
        couplingTerm

        helpers
        pmin   = 1e-1*barsa % cutoff value used in updateState method
        mfmin  = 0.01       % cutoff value used in updateState method
        mfmax  = 0.9        % cutoff value used in updateState method
        phimax = 5          % cutoff value used in updateState method
        
    end
    
    methods
        
        function model = ProtonicMembraneCell(inputparams)

            model = model@BaseModel();

            fdnames = {'G', ...
                       'T', ...
                       'couplingTerm'};
            
            model = dispatchParams(model, inputparams, fdnames);

            model.Electrolyser = ProtonicMembrane(inputparams.Electrolyser);
            model.GasSupply    = ProtonicMembraneGasSupply(inputparams.GasSupply);

            % Setup interface model
            interfaceinputparams.T                = inputparams.GasSupply.T;
            interfaceinputparams.molecularWeights = inputparams.GasSupply.molecularWeights;

            model.Interface = ProtonicMembraneGasSupplyBc(interfaceinputparams);
            
            model = model.setupInterface();
            
            model.constants = PhysicalConstants();

            model = model.setupScalings(inputparams);
            
        end
        
        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@BaseModel(model);

            gs     = 'GasSupply';
            elyser = 'Electrolyser';
            an     = 'Anode';
            ct     = 'Cathode';
            elyte  = 'Electrolyte';
            ctrl   = 'Control';
            itf    = 'Interface';
            
            nGas = model.GasSupply.nGas;

            varnames = {};

            % Coupling equations
            varnames{end + 1} = VarName({}, 'massCouplingEquations', nGas);
            % We add time variables to this model and the Electrolyser model. Needed for the control
            varnames{end + 1} = 'time';
            % Coupling damping variable
            varnames{end + 1} = 'beta';
            
            model = model.registerVarNames(varnames);

            % The partial pressures are not used directly in assembly but we want to compute them in post-processing
            model = model.unsetAsExtraVarName(VarName({'Interface'}, 'pressures', nGas));

            fn =  @ProtonicMembraneCell.updateInterfaceEquation;
            inputvarnames = {VarName({'Interface'}, 'massFluxes', nGas)   , ...
                             VarName({'Interface'}, 'massfractions', nGas), ...
                             VarName({'Interface'}, 'densities', nGas)    , ...
                             VarName({'Interface'}, 'pressure')           , ...
                             VarName({'Interface'}, 'density')            , ...
                             VarName({'GasSupply'}, 'massfractions', nGas), ...
                             VarName({'GasSupply'}, 'densities', nGas)    , ...
                             VarName({'GasSupply'}, 'pressure')           , ...
                             VarName({'GasSupply'}, 'density')};
            outputvarname = VarName({'Interface'}, 'bcFluxEquations', nGas);
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});

            fn =  @ProtonicMembraneCell.updateAnodePressures;
            inputvarnames = {VarName({'Interface'}, 'pressures', nGas)};
            outputvarname = VarName({'Electrolyser', 'Anode'}, 'pressures', nGas);
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});

            fn =  @ProtonicMembraneCell.updateCouplingEquation;
            inputvarnames = {{'Electrolyser', 'Anode' 'iHp'}, ...
                             VarName({'Interface'}, 'massFluxes', nGas), ...
                             'beta'};
            outputvarname = VarName({}, 'massCouplingEquations', nGas);
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});
            
            fn =  @ProtonicMembraneCell.updateGasSupplySource;
            inputvarnames = {VarName({'Interface'}, 'massFluxes', nGas), VarName({'GasSupply', 'GasSupplyBc'}, 'massFluxes', nGas)};
            outputvarname = VarName({'GasSupply'}, 'massSources', nGas);
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});
            
            model = model.setAsStaticVarName('time');
            
            fn =  @ProtonicMembraneCell.dispatchTime;
            model = model.registerPropFunction({{elyser, 'time'}, fn, {'time'}});

            fn = @ProtonicMembraneCell.updateBeta;
            fn = {fn, @(propfunction) PropFunction.drivingForceFuncCallSetupFn(propfunction)};
            model = model.registerPropFunction({'beta', fn, {'time'}});
            
        end


        function model = setupScalings(model, inputparams)

            jsonstruct = inputparams.jsonstruct;

            elyser = 'Electrolyser';
            elyte  = 'Electrolyte';
            an     = 'Anode';
            ct     = 'Cathode';
            ctrl   = 'Control';
            
            elysermodel = ProtonicMembrane(inputparams.(elyser));
            elysermodel = elysermodel.setupForSimulation();

            initElyserState = elysermodel.setupInitialState();

            initElyserState.Control.I = 0;

            drivingforces.src = @(time, I) pmControlFunc(time, I, 1, 2, 'order', 'I-first');
            initElyserState = elysermodel.evalVarName(initElyserState, {elyte, 'sigmaEl'}, {{'drivingForces', drivingforces}});
            initElyserState = elysermodel.evalVarName(initElyserState, {elyte, 'sigmaHp'}, {{'drivingForces', drivingforces}});

            sigmaHp = initElyserState.(elyte).sigmaHp(1);
            sigmaEl = initElyserState.(elyte).sigmaEl(1);
            
            phi0    = abs(model.(elyser).(an).E_0 - model.(elyser).(ct).E_0); % characteristic voltage
            T       = model.(elyser).(elyte).G.getTrans();
            T       = T(1);

            sHp = T*sigmaHp*phi0;
            sEl = T*sigmaEl*phi0;

            elyserScalings = {{{elyser, elyte, 'massConsHp'}     , sHp}, ...
                              {{elyser, elyte, 'chargeConsEl'}   , sEl}, ...
                              {{elyser, an   , 'chargeCons'}     , sEl}, ...
                              {{elyser, an   , 'iElEquation'}    , sEl}, ...
                              {{elyser, an   , 'iHpEquation'}    , sHp}, ...
                              {{elyser, ct   , 'chargeCons'}     , sEl}, ...
                              {{elyser, ct   , 'iElEquation'}    , sEl}, ...
                              {{elyser, ct   , 'iHpEquation'}    , sHp}, ...
                              {{elyser, ctrl , 'controlEquation'}, sEl}, ...
                              {{elyser, 'anodeChargeCons'}       , sEl}};

            schedule = elysermodel.Control.setupSchedule(jsonstruct.(elyser));

            [~, states, report] = simulateScheduleAD(initElyserState, elysermodel, schedule);

            state = states{end};
            state = elysermodel.addVariables(state, schedule.control);

            molfluxref = abs(sum(state.Anode.iHp)/PhysicalConstants.F); % in mol/second

            gs = 'GasSupply';
            
            gasInd = model.(gs).gasInd;
            pH2O = jsonstruct.(gs).initialState.pressure;


            switch jsonstruct.Geometry.type
              case '2D'
                scalFlux = molfluxref/model.(gs).molecularWeights(1)/jsonstruct.Geometry.Ny; % in kg/s
              otherwise
                error('geometry not covered yet');
            end
            
            scalPressure = pH2O;
            scalRateEq = jsonstruct.(gs).Nx*scalFlux;
            
            gsScalings = {{{gs, 'massConses', 1}, scalFlux}                    , ...
                          {{gs, 'massConses', 2}, scalFlux}                    , ...
                          {{gs, 'Control', 'pressureEq'}, scalPressure}        , ...
                          {{gs, 'Control', 'rateEq'}, scalRateEq}              , ...
                          {{gs, 'GasSupplyBc', 'bcFluxEquations', 1}, scalFlux}, ...
                          {{gs, 'GasSupplyBc', 'bcFluxEquations', 2}, scalFlux}};


            model.scalings = horzcat(elyserScalings, gsScalings);
            
        end
        
        function initstate = setupInitialState(model, jsonstruct)

            gs     = 'GasSupply';
            elyser = 'Electrolyser';
            an     = 'Anode';
            
            gasInd = model.(gs).gasInd;
            nGas   = model.(gs).nGas;
            Mws    = model.(gs).molecularWeights;
            
            model.(gs) = model.(gs).setupComputationalGraph();
            
            initstate.(gs) = model.(gs).setupInitialState(jsonstruct.GasSupply);

            p   = initstate.(gs).pressure;
            mfs = initstate.(gs).massfractions;
            
            tot = 0;
            for igas = 1 : nGas
                tot = tot + mfs{igas}/Mws(igas);
            end
            
            for igas = 1 : nGas
                pgs{igas} = p./tot.*mfs{igas}/Mws(igas);
            end

            pH2O = pgs{gasInd.H2O};
            pO2  = pgs{gasInd.O2};

            model.(elyser).(an).pH2O = pH2O(1);
            model.(elyser).(an).pO2  = pO2(1);
            
            model.(elyser) = model.(elyser).setupComputationalGraph();

            initstate.(elyser) = model.(elyser).setupInitialState();

            initstate = model.evalVarName(initstate, VarName({gs}, 'density'));
            initstate.time = 0;

            map = model.helpers.interfaceMap;

            initstate.Interface.pressure = map*p;
            initstate.Interface.massfractions{1} = map*mfs{1};
            
            for igas = 1 : nGas
               
                initstate.Interface.massFluxes{igas} = map*(0*pgs{igas});
                
            end
            
        end

        function newstate = addVariablesAfterConvergence(model, newstate, state)
            
            newstate = addVariablesAfterConvergence@BaseModel(model, newstate, state);

            % we add those values to be able to run addVariables method
            newstate.Electrolyser.Electrolyte.alpha = state.Electrolyser.Electrolyte.alpha;
            newstate.beta                           = state.beta;
            newstate.Electrolyser.Control.I         = state.Electrolyser.Control.I;
            
        end
            
        function state = addVariables(model, state)

        % Given a state where only the primary variables are defined, this
        % functions add all the additional variables that are computed in the assembly process and have some physical
        % interpretation.
        %
        % To do so, we use getEquations function and sends dummy variable for state0, dt and drivingForces

        % Values that need to be set to get the function getEquations running

            dt = 1;
            state0 = state;

            function [I, alpha, beta] = src(time, I)
            % we need only beta value
                I     = state.Electrolyser.Control.I;
                alpha = state.Electrolyser.Electrolyte.alpha;
                beta  = state.beta;
            end
            
            drivingForces.src = @(time, I) src(time, I);
            
            % We call getEquations to update state

            [~, state] = getEquations(model, state0, state, dt, drivingForces, 'ResOnly', true);

        end
        
        function model = setupInterface(model)

            coupterms = model.Electrolyser.couplingTerms;
            
            for icoup = 1 : numel(coupterms)
                coupterm = coupterms{icoup};
                [lia, locb] = ismember('Anode', coupterm.componentnames);
                if lia
                    loc = true(2, 1);
                    loc(locb) = false;
                    faces = coupterm.couplingfaces(:, loc);
                    break
                end
            end

            % face interface indices in parent grid.
            pfaces = model.Electrolyser.Electrolyte.G.mappings.facemap(faces);

            invfacemap = model.GasSupply.G.mappings.invfacemap;

            %  face interface indices in gas supply grid
            faces = invfacemap(pfaces);

            tbls = setupTables(model.GasSupply.grid);

            cellfacetbl = tbls.cellfacetbl;
            celltbl     = tbls.celltbl;
            
            bcfacetbl.faces = faces;
            bcfacetbl = IndexArray(bcfacetbl);
            bcfaceindtbl = bcfacetbl.addInd('ind', (1 : bcfacetbl.num)');

            bccellfaceindtbl = crossIndexArray(bcfaceindtbl, cellfacetbl, {'faces'});
            bccellfaceindtbl = sortIndexArray(bccellfaceindtbl, {'ind', 'faces', 'cells'});

            map = TensorMap();
            map.fromTbl = celltbl;
            map.toTbl = bccellfaceindtbl;
            map.mergefds = {'cells'};

            M = SparseTensor();
            M = M.setFromTensorMap(map);
            M = M.getMatrix();

            
            T = model.GasSupply.G.getBcTrans(bccellfaceindtbl.get('faces'));
            
            model.helpers.interfaceMap   = M;
            model.helpers.interfaceT     = T;
            model.helpers.interfaceCells = bccellfaceindtbl.get('cells');
            model.helpers.interfaceFaces = bccellfaceindtbl.get('faces');
            
        end


        function state = updateBeta(model, state, drivingForces)
            
            time = state.time;
            [~, ~, beta] = drivingForces.src(time, state.Electrolyser.Control.I);
            state.beta = beta;
            % state.beta = 0;
            
        end

        function state = dispatchTime(model, state)

            state.Electrolyser.time = state.time;
            
        end

        function state = updateInterfaceEquation(model, state)
            
        % Implementation is basically the same as updateBCequations method in ProtonicMembraneGasSupply, except for the
        % helper structure (hence an other implementation, but it could be cleaned-up and have just one implementation)
            
            map  = model.helpers.interfaceMap;
            Tbc  = model.helpers.interfaceT;
            nGas = model.Interface.nGas;
            
            K    = model.GasSupply.permeability;
            mu   = model.GasSupply.viscosity;
            D    = model.GasSupply.diffusionCoefficients;
            
            pIn = map*state.GasSupply.pressure;
            pBc = state.Interface.pressure;

            rhoBc = state.Interface.density;

            for igas = 1 : nGas

                mfBc      = state.Interface.massfractions{igas};
                bcFlux    = state.Interface.massFluxes{igas};
                rhoInigas = map*state.GasSupply.densities{igas};
                rhoBcigas = state.Interface.densities{igas};

                bceqs{igas} = rhoBc.*mfBc*K/mu.*Tbc.*(pIn - pBc) + D(igas).*Tbc.*(rhoInigas - rhoBcigas) - bcFlux;

            end

            state.Interface.bcFluxEquations = bceqs;

        end


        function state = updateGasSupplySource(model, state)

            nGas     = model.GasSupply.nGas;
            bccells  = model.GasSupply.helpers.bccells;
            intcells = model.helpers.interfaceCells;
            
            for igas = 1 : nGas

                src = 0.*state.GasSupply.pressure; % hacky initialization to get AD

                src(bccells)  = src(bccells) + state.GasSupply.GasSupplyBc.massFluxes{igas};
                src(intcells) = src(intcells) + state.Interface.massFluxes{igas};

                srcs{igas} = - src;
                
            end

            state.GasSupply.massSources = srcs;
            
        end

        
        function state = updateCouplingEquation(model, state)

            nGas   = model.GasSupply.GasSupplyBc.nGas;
            gasInd = model.GasSupply.GasSupplyBc.gasInd;
            mws    = model.GasSupply.molecularWeights;
            F      = PhysicalConstants.F;
            
            iHp     = state.Electrolyser.Anode.iHp;
            mfluxes = state.Interface.massFluxes;
            beta    = state.beta;
            
            % Chemical equation : 1/2*H2O <-> H^+ + e^- + 1/4*O2
            % Flux in the boundary conditions are oriented outwards from GasSupply domain.
            igas = gasInd.H2O;
            coupeqs{igas} = mfluxes{igas} - beta*1/2*mws(igas)*iHp/F;

            igas = gasInd.O2;
            coupeqs{igas} = mfluxes{igas} + beta*1/4*mws(igas)*iHp/F;

            state.massCouplingEquations = coupeqs;
            
        end

        function state = updateAnodePressures(model, state)

            state.Electrolyser.Anode.pressures = state.Interface.pressures;
            
        end


        function model = setupForSimulation(model)

            model.isRootSimulationModel = true;

            shortNames = {'1'                    , 'H2O';
                          '2'                    , 'O2';
                          'massConses'           , 'massCons';
                          'Electrolyser'         , 'elyser';
                          'GasSupply'            , 'gs';
                          'GasSupplyBc'          , 'gsBc';
                          'massCouplingEquations', 'massCoupEqs';
                          'controlEquations'     , 'ctrlEqs';
                          'bcFluxEquations'      , 'bcEqs';
                          'Electrolyte'          , 'elyte';
                          'Anode'                , 'an';
                          'Cathode'              , 'ct';
                          'Control'              , 'ctrl'};

            model = model.equipModelForComputation('shortNames', shortNames);

        end

        
        function forces = getValidDrivingForces(model)

            forces = getValidDrivingForces@BaseModel(model);
            forces.src   = [];
            
        end

        function model = validateModel(model, varargin)
        % do nothing
        end

        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            
            [state, report] = updateState@BaseModel(model, state, problem, dx, drivingForces);

            % cutoff the pressures

            gs     = 'GasSupply';
            elyser = 'Electrolyser';
            an     = 'Anode';
            ct     = 'Cathode';
            itf    = 'Interface';
            elyte  = 'Electrolyte';
            ctrl   = 'Control';

            pmin   = model.pmin;
            mfmin  = model.mfmin;
            mfmax  = model.mfmax;
            phimax = model.phimax;
            

            pressurevars = {{gs, 'pressure'}, ...
                            {gs, 'GasSupplyBc', 'pressure'}, ...
                            {itf, 'pressure'}};

            for ip = 1 : numel(pressurevars)
                state = model.capProperty(state, pressurevars{ip}, pmin);
            end

            mfvars = {{gs, 'massfractions', 1}               , ...
                      {gs, 'GasSupplyBc', 'massfractions', 1}, ...
                      {itf, 'massfractions', 1}};
            
            for ip = 1 : numel(mfvars)
                state = model.capProperty(state, mfvars{ip}, mfmin, mfmax);
            end
            
            phivars = {{'Electrolyser', 'Anode'      , 'phi'}, ...
                       {'Electrolyser', 'Cathode'    , 'pi'} , ...
                       {'Electrolyser', 'Electrolyte', 'phi'}, ...
                       {'Electrolyser', 'Electrolyte', 'pi'}};
            
            for ip = 1 : numel(phivars)
                state = model.capProperty(state, phivars{ip}, -phimax , phimax);
            end

        end


        function step = setupScheduleStep(model, jsonstruct)

            jtstruct = jsonstruct.TimeStepping;
            
            N = jtstruct.numberOfTimeSteps;
            T = jtstruct.totalTime;
            f = jtstruct.fractionSwitch;

            timeswitch = f*T;

            % we dedicate about 1/3 of time steps to the first phase
            N1 = floor(1/3*N);
            N2 = N - N1;
            
            dt1 = timeswitch/N1;
            t1  = linspace(0, dt1, N1 + 1)';
            alpha = 40; % we invert the exponential growth (to try to match with the exponential function in the electronic conductivity)
            t1 = log(alpha*t1 + 1)./log(alpha*dt1 + 1)*dt1;

            steps1 = diff(t1);

            dt2 = (T - timeswitch)/N2;
            steps2 = rampupTimesteps(totaltime - timeswitch, dt2, 5);

            step.val = [steps1; steps2];
            step.control = ones(numel(step.val), 1);

        end

        function control = setupScheduleControl(model, jsonstruct)

            jtstruct = jsonstruct.TimeStepping;

            T           = jtstruct.totalTime;
            f           = jtstruct.fractionSwitch;
            orderSwitch = jtstruct.orderSwitch;

            timeswitch = f*T;
            
            control.src = @(time, I) pmControlFunc(time, I, timeswitch, T, 'order', orderSwitch);
            
        end
        
        function schedule = setupSchedule(model, jsonstruct)
        % Convenience function to setup schedule from main jsonstruct with property TimeStepping

            step    = model.setupScheduleStep(jsonstruct);
            control = model.setupScheduleControl(jsonstruct);
            
            schedule = struct('step', step, ...
                              'control', control);
            
        end
        
        
    end
    
end
