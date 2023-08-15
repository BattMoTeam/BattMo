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
