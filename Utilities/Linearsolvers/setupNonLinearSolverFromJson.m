function [model, nls, jsonstruct] = setupNonLinearSolverFromJson(model, jsonstruct)

    %% Setup the properties of the nonlinear solver 
    nls = NonLinearSolver();

    % setup default values
    jsonstruct = setDefaultJsonStructField(jsonstruct, {'NonLinearSolver', 'maxIterations'},  10);
    jsonstruct = setDefaultJsonStructField(jsonstruct, {'NonLinearSolver', 'maxTimestepCuts'},  6);
    jsonstruct = setDefaultJsonStructField(jsonstruct, {'NonLinearSolver', 'nonlinearTolerance'},  []);
    jsonstruct = setDefaultJsonStructField(jsonstruct, {'NonLinearSolver', 'verbose'},  false);
    jsonstruct = setDefaultJsonStructField(jsonstruct, {'NonLinearSolver', 'timeStepSelector'}, 'stateChangeTimeStepSelector');

    linearSolverSetup_default.library = 'matlab';
    linearSolverSetup_default.method  = 'direct';
    jsonstruct = setDefaultJsonStructField(jsonstruct, {'NonLinearSolver', 'LinearSolver', 'linearSolverSetup'}, linearSolverSetup_default);
    linearSolverSetup = getJsonStructField(jsonstruct, {'NonLinearSolver', 'LinearSolver', 'linearSolverSetup'});

    if ~strcmp(linearSolverSetup.library, 'matlab') || ~strcmp(linearSolverSetup.method, 'direct')
        assert(BatteryLinearSolver.checkAMGCL(), 'AMGCL not installed correctly');
    end
    
    nls.maxIterations   = jsonstruct.NonLinearSolver.maxIterations;
    nls.maxTimestepCuts = jsonstruct.NonLinearSolver.maxTimestepCuts;
    nls.verbose         = jsonstruct.NonLinearSolver.verbose;
    
    % Change default behavior of nonlinear solver, in case of error
    nls.errorOnFailure    = false;
    nls.continueOnFailure = false;

    switch getJsonStructField(jsonstruct, {'NonLinearSolver', 'timeStepSelector'})
      case 'stateChangeTimeStepSelector'
        nls.timeStepSelector = StateChangeTimeStepSelector('TargetProps', {{'Control','E'}}, 'targetChangeAbs', 0.03);
      case 'simpleTimeStepSelector'
        nls.timeStepSelector = SimpleTimeStepSelector();
      otherwise
       error('Unknown timeStepSelector');
    end
    
    nls.LinearSolver = BatteryLinearSolver('linearSolverSetup', linearSolverSetup);
    
    % Change default tolerance for nonlinear solver
    % For the moment, this is hard-coded here
    if isprop(model.Control, 'Imax')
        ImaxRef = model.Control.Imax;
    elseif isprop(model.Control, 'ImaxDischarge')
        ImaxRef = model.Control.ImaxDischarge;
    else
        ImaxRef = [];
    end
        
    if isfield(jsonstruct.NonLinearSolver, 'nonlinearTolerance') && ~isempty(jsonstruct.NonLinearSolver.nonlinearTolerance)
        model.nonlinearTolerance = jsonstruct.NonLinearSolver.nonlinearTolerance;
    elseif ~isempty(ImaxRef)
        model.nonlinearTolerance = 1e-3*ImaxRef;
    else
        % some default value
        model.nonlinearTolerance = 1e-5;
    end
    
    % Set verbosity
    model.verbose = jsonstruct.NonLinearSolver.verbose;

    if isfield(linearSolverSetup, 'reduction') && linearSolverSetup.reduction.doReduction
        model = model.setupSelectedModel('reduction', linearSolverSetup.reduction);
    end
        
end



%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
