classdef seiNonLinearSolver < NonLinearSolver

    methods
        function solver = seiNonLinearSolver(varargin)
            solver = solver@NonLinearSolver(varargin{:});
            
        end
        
        function [state, failure, report] = solveMinistep(solver, model, state, state0, dt, drivingForces)
            
            [state, failure, report] = solveMinistep@NonLinearSolver(solver, model, state, state0, dt, drivingForces);
            
            if report.Converged
                if strcmp(state.Control.ctrlType, 'CC_discharge1')
                    Emin = model.Control.lowerCutoffVoltage;
                    E = state.Control.E;
                    if  E < Emin
                        failure = true;
                        report.Converged = false;
                    end
                end
            end
            
        end
        
    end
    
end
