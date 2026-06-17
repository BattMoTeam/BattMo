classdef ElectrolyserImpedanceSolver < ImpedanceSolver

    methods

        function impsolv = ElectrolyserImpedanceSolver(inputparams, options, extrastructs)

            if nargin < 3
                extrastructs = [];
                if nargin < 2
                    options = [];
                end
            end

            impsolv = impsolv@ImpedanceSolver(inputparams, options, extrastructs);
            
            % other choice for state initializaztion is 'given state'
            options = setDefaultJsonStructField(options, {'stateInitialization', 'initializationSetup'}, 'given current');

            if strcmp(getJsonStructField(options, {'stateInitialization', 'initializationSetup'}), 'given current')
                options = setDefaultJsonStructField(options, {'stateInitialization', 'computeSteadyState'}, true);
            end
            
        end

        function setupIvalue(impsolv)

            options = impsolv.options;
            
            switch getJsonStructField(options, {'stateInitialization', 'initializationSetup'})
                
              case 'given current'

                impsolv.I = getJsonStructField(options, {'stateInitialization', 'I'});
                
              case 'given state'

                initstate = impsolv.extrastructs.initstate;
                impsolv.I = impsolv.model.getProp(initstate, impsolv.Ivarname);
                
              otherwise
                
                error('initializationSetup not recognized');
                
            end
            
        end
        
        function setupModel(impsolv)

            model = Electrolyser(impsolv.inputparams);
            model = model.setupBcAndInitialState();

            impsolv.model = model;
            
        end
        
        function setupVarNames(impsolv)
            
            oer = 'OxygenEvolutionElectrode';
            ptl = 'PorousTransportLayer';
            ctl = 'CatalystLayer';
            impsolv.Ivarname = {oer, ctl, 'I'};
            impsolv.Uvarname = {oer, ptl, 'E'};
            
        end


        function drivingForces = setupDrivingForces(impsolv)
            
            drivingForces = impsolv.model.getValidDrivingForces();
            drivingForces.src = @(time) impsolv.I;
            
        end

        function setupSteadyState(impsolv)

            options = impsolv.options;

            if getJsonStructField(options, {'stateInitialization', 'computeSteadyState'})

                model    = impsolv.model;
                controlI = impsolv.I;
                
                switch getJsonStructField(options, {'stateInitialization', 'initializationSetup'})
                    
                  case 'given current'

                    [~, initstate] = model.setupBcAndInitialState();
                    
                  case 'given state'

                    initstate = impsolv.extrastructs.initstate;
                    
                  otherwise
                    
                    error('initializationSetup not recognized');
                    
                end
                
                total = 10*hour;

                n   = 10;
                dt  = total/n;
                dts = rampupTimesteps(total, dt, 5);

                tup      = total;
                srcfunc  = @(time) rampupControl(time, tup, controlI, 'rampupcase', 'linear');
                control  = struct('src', srcfunc);

                step = struct('val', dts, 'control', ones(numel(dts), 1));
                schedule = struct('control', control, 'step', step);

                nls = NonLinearSolver();
                
                nls.verbose        = false;
                nls.errorOnFailure = false;
                
                simsetupInput = struct('model'          , model    , ...
                                       'initstate'      , initstate, ...
                                       'NonLinearSolver', nls      , ...
                                       'schedule'       , schedule);

                %% Run the simulation
                
                simsetup.model.verbose = true;

                simsetup = SimulationSetup(simsetupInput);
                
                states = simsetup.run();

                impsolv.state = states{end};

            else
                
                impsolv.state = impsolv.extrastructs.initstate;
                
            end
        end
        
    end
    
end

