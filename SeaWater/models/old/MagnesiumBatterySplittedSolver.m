classdef MagnesiumBatterySplittedSolver < BaseModel

    properties
        
        noPrecModel
        noFlowModel
        fullModel
        
    end
    
    methods

        function model = MagnesiumBatterySplittedSolver(paramobj)
            
            model = model@BaseModel();
            
            model.noPrecModel = MagnesiumBatteryNoPrecipitation(paramobj);
            model.noFlowModel = MagnesiumBatteryNoFlow(paramobj);
            model.fullModel   = MagnesiumBattery(paramobj);

        end

        
        function [state, report] = stepFunction(model, state, state0, dt, drivingForces, linsolver, nonlinsolver, iteration, varargin)
            
            if iteration == 0
                
                for ind = 1
                    
                    donoflowiteration = true;
                    
                    if donoflowiteration
                        
                        fprintf('\n** Starting NO FLOW solver step - iteration %d.\n\n', iteration);
                        [state, failure, noflowreport] = nonlinsolver.solveMinistep(model.noFlowModel, state, state0, dt, drivingForces);
                        
                        if ~noflowreport.Converged
                            report = model.makeStepReport();
                            report.Failure = true;
                            report.FailureMsg = 'no flow solver did not converge';
                            report.Converged = false;
                            return
                        end
                        
                    end

                    dopureflowiteration = true;
                    
                    if dopureflowiteration
                        fprintf('\n** Starting NO PRECIPITATION solver step.\n\n')

                        [state, failure, noprecreport] = nonlinsolver.solveMinistep(model.noPrecModel, state, state0, dt, drivingForces);
                        
                        if ~noprecreport.Converged
                            report = model.makeStepReport();
                            report.Failure = true;
                            report.FailureMsg = 'no prec solver did not converge';
                            report.Converged = false;
                            fprintf('\n** NO PRECIPITATION solver step did not converge.\n\n')
                            return
                        end                

                    end
                    
                end
                
                fprintf('\n** Starting FULL MODEL step.\n\n')
                
            end
            
            [state, report] = model.fullModel.stepFunction(state, state0, dt, drivingForces, linsolver, nonlinsolver, iteration, varargin{:});
            
        end

        
        function validforces = getValidDrivingForces(model)
            validforces=struct('src', [], 'stopFunction', []);
        end
        
        function model = validateModel(model, varargin)
            model.fullModel = model.fullModel.validateModel(varargin{:});
            model.noPrecModel = model.noPrecModel.validateModel(varargin{:});
        end
        
        function cleanState = addStaticVariables(model, cleanState, state)
            cleanState = addStaticVariables@BaseModel(model, cleanState, state);
            cleanState.time = state.time;
        end
        
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)

            [state, report] = model.fullModel.updateAfterConvergence(state0, state, dt, drivingForces);
            
        end

                
    end

end

