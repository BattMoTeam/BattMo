function [model, nls] = setupNonLinearSolverFromJson(model, jsonstruct)

    %% Setup the properties of the nonlinear solver 
    nls = NonLinearSolver();

    % load from json file or use default value. 
    if ~isfield(jsonstruct, 'NonLinearSolver')
        
        % setup default values
        jsonstruct.NonLinearSolver.maxIterations = 10;
        jsonstruct.NonLinearSolver.nonlinearTolerance = [];
        jsonstruct.NonLinearSolver.verbose = false;
        
        linearSolverSetup.library = "matlab";
        linearSolverSetup.method = "direct";

    else

        linearSolverSetup = jsonstruct.NonLinearSolver.LinearSolver.linearSolverSetup;
        
    end
    
    nls.maxIterations = jsonstruct.NonLinearSolver.maxIterations;
    if isfield(jsonstruct.NonLinearSolver, 'verbose')
        nls.verbose = jsonstruct.NonLinearSolver.verbose;
    else
        nls.verbose = false;
    end
    
    % Change default behavior of nonlinear solver, in case of error
    nls.errorOnFailure = false;
    nls.timeStepSelector = StateChangeTimeStepSelector('TargetProps', {{'Control','E'}}, 'targetChangeAbs', 0.03);

    nls.LinearSolver = BatteryLinearSolver('linearSolverSetup', linearSolverSetup);
    
    % Change default tolerance for nonlinear solver
    % For the moment, this is hard-coded here
    if isfield(jsonstruct.NonLinearSolver, 'nonlinearTolerance') && ~isempty(jsonstruct.NonLinearSolver.nonlinearTolerance)
        model.nonlinearTolerance = jsonstruct.NonLinearSolver.nonlinearTolerance;
    elseif ~isempty(model.Control.Imax)
        model.nonlinearTolerance = 1e-3*model.Control.Imax;
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
