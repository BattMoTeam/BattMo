classdef BatteryLinearSolver < handle
% Most of the basic methods are essentially taken from MRST LinearSolver class
    
    properties

        tolerance                  % Linear solver tolerance
        maxIterations              % Max number of iterations used
        extraReport                % Enable this to produce additional report output
                                   % May use a lot of memory for large problems
        verbose                    % Verbose output enabler
        replaceNaN                 % Boolean indicating if the solver should replace NaN in the results
        replaceInf                 % Boolean indicating if the solver should replace Inf in the results
        replacementNaN             % If replaceNaN is enabled, this is the value that will be inserted
        replacementInf             % If replaceInf is enabled, this is the value that will be inserted
        reduceToCell               % Reduce to per-cell system before solving
        applyLeftDiagonalScaling   % Apply left diagonal scaling before solving
        applyRightDiagonalScaling  % Apply right diagonal scaling before solving
        keepNumber                 % If set, linear solver will reduce the system to the first keepNumber entries
        useSparseReduction         % If true, sparse indexing will be used with keepNumber option
        variableOrdering           % Variable ordering to be used for linear solver. Row vector of equal length to the size of the linear system.
        equationOrdering           % Equation ordering to be used for linear solver. Row vector of equal length to the size of the linear system.
        id = '';                   % Short text string identifying the specific solver. Appended to the short name (see getDescription)
        initialGuessFun = [];

        setup  % Structure described in schema battmodDir()/Utilities/JsonSchemas/linearsolver.schema.json (recommended reference)
               %
               % the first fields are
               % - method : string which can be any of
               %            - 'direct'    : no iterative solver, we use matlab direct inversion
               %            - 'iterative' : We use iterative solver for the whole system without preconditioning.
               %                            The solver description is sent into structure solver (see examples below)
               %            - 'gmres'     : We use matlab gmres solver with preconditioning. Further setup is sent into options structure (see example below)
               % QUESTION : add more here
        
        verbosity
        first
        reuse_setup
        precondIterations_phi
        precondIterations_c
        precondIterations_T
        
    end

    methods
        
        function solver = BatteryLinearSolver(varargin)
            
            solver.tolerance                 = 1e-8;
            solver.maxIterations             = 25;
            solver.extraReport               = false;
            solver.verbose                   = mrstVerbose();
            solver.replaceNaN                = false;
            solver.replaceInf                = false;
            solver.replacementNaN            = 0;
            solver.replacementInf            = 0;
            solver.reduceToCell              = false;
            solver.useSparseReduction        = false;
            solver.applyLeftDiagonalScaling  = false;
            solver.applyRightDiagonalScaling = false;
            solver.variableOrdering          = [];
            solver.equationOrdering          = [];
            solver.verbosity                 = 0;
            solver.setup                     = struct('method', 'direct');
            solver.reuse_setup               = false;
            solver.first                     = true;
            
            solver = merge_options(solver, varargin{:});
            
            assert(solver.maxIterations >= 0);
            assert(solver.tolerance >= 0);
            
        end

        function [dx, result, report] = solveLinearProblem(solver, problem, model)
        % Solve a linearized problem
        % QUESTION : Do we want to simplify solveLinearProblem ?
            
            timer = tic();
            lsys = [];
            eliminated = {};
            keepNumber0 = solver.keepNumber;
            initialGuess = solver.getInitialGuess(problem);
            
            if solver.reduceToCell && isempty(solver.keepNumber)
                % Eliminate non-cell variables
                
                s = getSampleAD(problem.equations{:});
                keep = problem.indexOfType('cell');
                
                if ~all(keep)
                    
                    if isa(s, 'GenericAD')
                        % If we are working with block AD, we use the built-in
                        % keepNumber property of the linear solver to perform a
                        % full block Schur complement
                        nk = sum(keep);
                        assert(all(keep(1 : nk)) & ~any(keep(nk + 1 : end)), ...
                               'Cell variables must all combine first in the ordering for this AutodiffBackend.');
                        if 1
                            % In-place Schur complement
                            ngroups = numel(s.offsets)-1;
                            varno = rldecode((1:ngroups)', diff(s.offsets));

                            keepEq = problem.equations(keep);
                            elimEq = problem.equations(~keep);
                            keepVar = rldecode((keep)', diff(s.offsets));
                            keepVar = unique(varno(keepVar));
                            elimVar = setdiff(varno, keepVar);

                            B_eq = keepEq;
                            C_eq = B_eq;
                            for i = 1:numel(B_eq)
                                B_eq{i}.jac = B_eq{i}.jac(keepVar);
                                C_eq{i}.jac = C_eq{i}.jac(elimVar);
                            end
                            D_eq = elimEq;
                            E_eq = D_eq;
                            for i = 1:numel(D_eq)
                                D_eq{i}.jac = D_eq{i}.jac(keepVar);
                                E_eq{i}.jac = E_eq{i}.jac(elimVar);
                            end
                            B_eq = combineEquations(B_eq{:});
                            C_eq = combineEquations(C_eq{:});
                            D_eq = combineEquations(D_eq{:});
                            E_eq = combineEquations(E_eq{:});
                            lsys = struct('B', B_eq.jac{1}, ...
                                          'C', C_eq.jac{1}, ...
                                          'D', D_eq.jac{1}, ...
                                          'E', E_eq.jac{1},...
                                          'f', -B_eq.val, ...
                                          'h', -D_eq.val, ...
                                          'E_L', [], ...
                                          'E_U', []);
                            [lsys.E_L, lsys.E_U] = lu(lsys.E);
                            problem.A = lsys.B - lsys.C*(lsys.E_U\(lsys.E_L\lsys.D));
                            problem.b = lsys.f - lsys.C*(lsys.E_U\(lsys.E_L\lsys.h));
                        else
                            nv =  s.getNumVars();
                            solver.keepNumber = sum(nv(keep));
                        end
                        
                    else
                        
                        [problem, eliminated] = solver.reduceToVariable(problem, keep);
                        
                    end
                    
                    if ~isempty(initialGuess{1})
                        initialGuess = initialGuess(keep);
                    end
                end
            end
            
            problem = problem.assembleSystem();
            if(not(all(isfinite(problem.b))))
                assert(all(isfinite(problem.b)), 'Linear system rhs must have finite entries.');
            end
            
            % Get linearized system
            [A, b] = problem.getLinearSystem();
            x0     = vertcat(initialGuess{:});

            % Reduce system (if not already done)
            if isempty(lsys)
                % QUESTION : why?
                [A, b, lsys, x0] = solver.reduceLinearSystem(A, b, false, x0);
            end
            % Reorder linear system
            [A, b, ordering, x0] = solver.reorderLinearSystem(A, b, [], x0);
            % Apply scaling
            [A, b, scaling, x0] = solver.applyScaling(A, b, x0);

            t_prepare = toc(timer);
            % Solve the system
            [result, report] = solver.solveLinearSystem(A, b, x0, problem);
            t_solve = toc(timer) - t_prepare;
            % Undo scaling
            result = solver.undoScaling(result, scaling);
            % Permute system back
            result = solver.deorderLinearSystem(result, ordering);
            % Recover eliminated variables on linear level
            result = solver.recoverLinearSystem(result, lsys);
            
            [result, report] = problem.processResultAfterSolve(result, report);
            report.SolverTime         = toc(timer);
            report.LinearSolutionTime = t_solve;
            report.PreparationTime    = t_prepare;
            report.PostProcessTime    = report.SolverTime - t_solve - t_prepare;
            
            if solver.replaceNaN
                result(isnan(result)) = solver.replacementNaN;
            end

            if solver.replaceInf
                result(isinf(result)) = solver.replacementInf;
            end

            dx = solver.storeIncrements(problem, result);

            if ~isempty(eliminated)
                dx = solver.recoverResult(dx, eliminated, keep);
            end
            
            solver.keepNumber = keepNumber0;
            
        end

        function [result, report] = solveLinearSystem(solver, A, b, x0, problem) % ok
            
            report = solver.getSolveReport();

            method = solver.setup.method;
            setup  = solver.setup;
            
            switch method
                
              case 'direct'                    

                result = A\b;

              case 'iterative'

                % QUESTION : We could add the solver option in method.options
                % QUESTION : We could add "amgcl"

                % QUESTION : Agree on assert?
                assert(solver.reduceToCell, 'agmg makes only sense if cell reduction has been done');
                assert(strcmp(setup.solverspec.name, 'agmg'), 'amgcl not supported here');
                
                a = tic();
                
                if (solver.reuse_setup)

                    if (solver.first)
                        % QUESTION : what does solver.first mean?
                        
                        solver.first = false;
                        agmg(A, b, 20, solver.tolerance, solver.maxIterations, solver.verbosity, [], -1); 
                        agmg(A, b, 20, solver.tolerance, solver.maxIterations, solver.verbosity, [], 1);
                        
                    end
                    
                    [result, flag, relres, iter] = agmg(A, b, 20, solver.tolerance, solver.maxIterations, solver.verbosity, [], 2);
                    
                    if (flag == 1)
                        
                        agmg(A, b, 20, solver.tolerance, solver.maxIterations, solver.verbosity, [], -1); 
                        agmg(A, b, 20, solver.tolerance, solver.maxIterations, solver.verbosity, [], 1); 
                        [result, flag, relres, iter_new] = agmg(A, b, 20, solver.tolerance, solver.maxIterations, solver.verbosity, result, 2); 
                        iter = iter + iter_new; 
                        solver.first = true;
                        
                    end
                    
                else
                    
                    [result, flag, relres, iter] = agmg(A, b, 20, solver.tolerance, solver.maxIterations, solver.verbosity);
                    
                end
                
                report.Iterations         = iter;
                report.Residual           = relres;
                report.Converged          = flag;
                report.LinearSolutionTime = toc(a);
                
              case 'gmres'

                assert(solver.reduceToCell, 'agmg makes only sense if cell reduction has been done');

                options = setup.options;

                switch options.method
                    
                  case 'grouped'

                    % QUESTION : reset works in this case ?
                    agmg(A, b, 20, solver.tolerance, solver.maxIterations, solver.verbosity, [], -1); 
                    agmg(A, b, 20, solver.tolerance, solver.maxIterations, solver.verbosity, [], 1);
                    
                    % QUESTION : error in argument?
                    % QUESTION : again we can send option there
                    f = @(b) solver.agmgprecond(b, A)
                    % f = @(b) solver.agmgprecond(A, b)
                    
                    a=tic;

                    [result, flags, relres, iter] = gmres(A, b, solver.maxIterations, solver.tolerance, solver.maxIterations, f);

                    report.Iterations         = (iter(1) - 1)*solver.maxIterations + iter(2);
                    report.Residual           = relres;
                    report.Converged          = flags;
                    report.LinearSolutionTime = toc(a);
                    
                  case 'separate'

                    solvers = options.solvers;
                    
                    solver.precondIterations_phi = 0;
                    solver.precondIterations_c   = 0;
                    solver.precondIterations_T   = 0;
                    
                    indphi = solver.getPotentialIndex(problem);
                    indc   = solver.getCIndex(problem);
                    indT   = solver.getTIndex(problem);
                    
                    lverbosity = 0;
                    ltol = 1e-3;
                    lmax = 40;

                    % QUESTION : We could send separate options for the solver
                    for isolver = 1 : numel(solvers)
                        varsolver = solvers{isolver};
                        switch varsolver.variable
                          case 'phi'
                            phi_solver = solver.getElipticSolver(varsolver.solverspec);
                          case 'c'
                            c_solver   = solver.getElipticSolver(varsolver.solverspec);
                          case 'T'
                            T_solver   = solver.getElipticSolver(varsolver.solverspec);
                        end
                    end
                    
                    assert(all(indphi + indc + indT == 1), 'We have some unknown in the system that are not expected') 

                    if ~any(indT)
                        % handle the case with no temperature
                        T_solver = [];
                    end
                    
                    precond = @(b) solver.p_gs_precond(b, A, indphi, indc, indT, phi_solver, c_solver, T_solver, []);
                    
                    a=tic;
                    restart = solver.maxIterations;

                    [result, flags, relres, iter] = gmres(A, b, restart, solver.tolerance, solver.maxIterations, precond);
                    
                    % add diagnostic fields in report
                    report.Iterations            = (iter(1) - 1)*solver.maxIterations + iter(2);
                    report.Residual              = relres;
                    report.Converged             = flags;
                    report.precondIterations_phi = solver.precondIterations_phi; 
                    report.precondIterations_c   = solver.precondIterations_c;
                    report.precondIterations_T   = solver.precondIterations_T;
                    report.LinearSolutionTime    = toc(a);

                  otherwise

                    error('choice not recognized')
                end

                
              case 'matlab_cpr_agmg'

                % QUESTION : do we keep that case? If yes, fix below
                
                indb = solver.getPotentialIndex(problem);
                % QUESTION : is agmg reset needed here?
                agmg(A(indb, indb), b(indb), 20, solver.tolerance, solver.maxIterations, solver.verbosity, [], -1); 
                agmg(A(indb, indb), b(indb), 20, solver.tolerance, solver.maxIterations, solver.verbosity, [], 1);

                %rphi=  agmg(A(indb,indb),b(indb),20,1e-4,20,1,,-1);
                %rphi=  agmg(A(indb,indb),b(indb),20,1e-4,20,1,[],1);
                solver_type = 'agmgsolver';
                error('next call is now not supported')
                voltage_solver = solver.getElipticSolver(solver_type, A, indb)
                
                inds=true(size(indb));%not(indb);
                                      %inds = not(indb);
                smoother_type = 'gs'                      
                smoother = getSmoother(smoother_type);                     
                f=@(b) solver.precondcpr(b,A,indb,voltage_solver,smoother)
                tic;
                result = gmres(A,b,25,solver.tolerance,solver.maxIterations,f);
                toc;
                
              case 'amgcl'

                % QUESTION : do we keep that case ?
                
                precondcase = 'ILU0';
                
                switch precondcase
                    
                  case 'ILU0'
                    
                    isolver = struct('type'     ,'gmres', ...
                                     'M'        ,200    , ...
                                     'verbosity',10);

                    precond = struct('class'    ,'relaxation', ...
                                     'type'     ,'ilu0'      , ...
                                     'damping'  ,1           , ...
                                     'verbosity',10);

                    options = struct('solver'      ,isolver  , ...
                                     'precond'     ,precond  , ...
                                     'solver_type' ,'regular', ...
                                     'write_params',true     , ...
                                     'block_size'  ,1        , ...
                                     'verbosity'   ,0        , ...
                                     'reuse_mode'  ,1);
                  case 'amg'
                    
                    %isolver = struct('type','gmres','M',20,'verbosity',3);
                    isolver = struct('type', 'bicgstab', ...
                                     'M'   , 50);
                    relaxation=struct('type', 'ilu0');
                    %relaxation=struct('type','spai0')
                    %relaxation=struct('type',lower(smoother))
                    %maxNumCompThreads(np)
                    alpha=0.001;
                    aggr=struct('eps_strong', alpha);
                    % coarsening=struct('type','aggregation','over_interp',1.0,'aggr',aggr)
                    %coarsening=struct('type','smoothed_aggregation','relax',2/3,'aggr',aggr,'estimate_spectral_radius',false,'power_iters',10,'over_interp',1.0)
                    coarsetarget=1200;

                    coarsening=struct('type'         , 'ruge_stuben', ...
                                      'rs_eps_strong', alpha        , ...
                                      'rs_trunc'     , true         , ...
                                      'rs_eps_trunc' , alpha);

                    precond = struct('class'        , 'amg'       , ...
                                     'coarsening'   , coarsening  , ...
                                     'relax'        , relaxation  , ...
                                     'coarse_enough', coarsetarget, ...
                                     'max_levels'   , 20          , ...
                                     'ncycle'       , 1           , ...
                                     'npre'         , 1           , ...
                                     'npost'        , 1           , ...
                                     'pre_cycle'    , 0           , ...
                                     'direct_coarse', true);
                    
                    options = struct('solver'      , isolver  , ...
                                     'precond'     , precond  , ...
                                     'reuse_mode'  , 1        , ...
                                     'solver_type' , 'regular', ...
                                     'write_params', false    , ...
                                     'block_size'  , 1        , ...
                                     'verbosity'   , 10);
                    
                  case 'amg_cpr'
                    
                    %isolver = struct('type','gmres','M',20,'verbosity',3);
                    isolver = struct('type', 'bicgstab', ...
                                     'M'   , 50);
                    relaxation=struct('type','ilu0')
                    %relaxation=struct('type','spai0')
                    %relaxation=struct('type',lower(smoother))
                    %maxNumCompThreads(np)
                    alpha=0.001;
                    aggr=struct('eps_strong', alpha);
                    % coarsening=struct('type','aggregation','over_interp',1.0,'aggr',aggr)
                    % coarsening=struct('type','smoothed_aggregation','relax',2/3,'aggr',aggr,'estimate_spectral_radius',false,'power_iters',10,'over_interp',1.0)
                    coarsetarget=1200;

                    coarsening=struct('type'         , 'ruge_stuben', ...
                                      'rs_eps_strong', alpha        , ...
                                      'rs_trunc'     , true         , ...
                                      'rs_eps_trunc' , alpha);
                    
                    precond = struct('class'        , 'amg'       , ...
                                     'coarsening'   , coarsening  , ...
                                     'relax'        , relaxation  , ...
                                     'coarse_enough', coarsetarget, ...
                                     'max_levels'   , 20          , ...
                                     'ncycle'       , 1           , ...
                                     'npre'         , 1           , ...
                                     'npost'        , 1           , ...
                                     'pre_cycle'    , 0           , ...
                                     'direct_coarse', true);
                    
                    options = struct('solver'      , isolver  , ...
                                     'precond'     , precond  , ...
                                     'reuse_mode'  , 1        , ...
                                     'solver_type' , 'regular', ...
                                     'write_params', false    , ...
                                     'block_size'  , 1        , ...
                                     'verbosity'   , 10);
                    
                  otherwise

                    error('precondcase not recognized');
                    
                end
                
                [result, extra] = amgcl(A, b, 'amgcloptions', options         , ...
                                              'blocksize'   , 1               , ...
                                              'tol'         , solver.tolerance, ...
                                              'maxiter'     , solver.maxIterations);
                
              otherwise
                
                error('Method not implemented');
                
            end
            
        end
        
        function indb = getPotentialIndex(solver, problem) % ok
        % QUESTION : here we see that diag structure is needed
        %QUESTION  : why include 'E' ?

            varinds1 = solver.getVarIndex(problem, {'*', 'phi'});
            varinds2 = solver.getVarIndex(problem, {'Control', 'E'});
            varinds = [varinds1, varinds2];
            indb = solver.getIndex(problem, varinds);
            
        end

        function indb = getCIndex(solver, problem) % ok

            varinds = solver.getVarIndex(problem, {'Electrolyte', 'c'});
            indb = solver.getIndex(problem, varinds);
            
        end

        function indb = getTIndex(solver, problem) % ok

            varinds = solver.getVarIndex(problem, {'ThermalModel', 'T'});
            indb = solver.getIndex(problem, varinds);
            
        end 

        function varinds = getVarIndex(solver, problem, varname)

            varinds = [];
            anymodel = strcmp('*', varname{1});
            
            for ivar = 1 : numel(problem.primaryVariables)
                pvarname = problem.primaryVariables{ivar};
                takeivar = true;

                if anymodel
                    if ~strcmp(varname{end}, pvarname{end})
                        takeivar = false;
                    end
                else
                    if numel(pvarname) == numel(varname)
                        iname = 1;
                        while takeivar & iname <= numel(varname)
                            if ~strcmp(pvarname{iname}, varname{iname})
                                takeivar = false;
                            end
                            iname = iname + 1;
                        end
                    else
                        takeivar = false;
                    end
                end
                if takeivar
                    varinds = [varinds, ivar];
                end

            end
        end
        
        function indb = getIndex(solver, problem, varinds)
            
            numVars = problem.equations{1}.getNumVars(); 
            % QUESTION : why call getNumVars for equations{1}?
            pos = cumsum(numVars); 
            pos = [[1; pos(1:end - 1) + 1], pos]; 
            posvar = pos(varinds, :); 
            ind = mcolon(posvar(:, 1), posvar(:, 2)); 
            indb = false(pos(end), 1); 
            indb(ind) = true;
            
        end
        
        function report = getSolveReport(solver, varargin) % ok
            report = struct('Iterations'        , 0, ... % Number of iterations (if iterative)
                            'Residual'          , 0, ... % Final residual
                            'SolverTime'        , 0, ... % Total time in solver
                            'LinearSolutionTime', 0, ... % Time spent solving system
                            'PreparationTime'   , 0, ... % Schur complement, scaling     , ...
                            'PostProcessTime'   , 0, ... % Recovery , undo scaling, ...
                            'Converged'         , true); % Bool indicating convergence
            report = merge_options_relaxed(report, varargin);
        end
        
        function r = p_gs_precond(solver, x, A, indphi, indc, indT, phi_solver, c_solver, T_solver, opt)% ok
        % QUESTION : find better name?
        % QUESTION : switch argument input x and A
        % QUESTION : support for without temperature?
            
            r = x*0; 
            % xp agmg(A(ind, ind), x(ind), 20, 1e-4, 20, 1, [], 2); 
            % lverb = false
            % xs = zeros(sum(not(indphi)), 1); 
            ldisp = @(tt) disp(tt); 
            % ldisp = @(tt) disp(''); 
            
            nr = 1; % number of GS iterations
            
            for kk = 1 : nr
                
                ldisp('voltage start')
                xp = x(indphi) - 0.0*A(indphi, not(indphi))*r(not(indphi));
                if(true)
                    [rp, flag, res, iter] = phi_solver(A(indphi, indphi), xp); 
                    solver.precondIterations_phi = solver.precondIterations_phi + iter; 
                else
                    ii = ind; 
                    rp = agmg(A(ii, ii), xp, 20, solver.tolerance, solver.maxIterations, solver.verbosity, [], 0); 
                end
                r(indphi) = rp; 
                ldisp('voltage end')

                xs = x(indc) - 0.0*A(indc, not(indc))*r(not(indc)); 

                ldisp('c start')
                if(true)
                    [rs, flag, res, iter] = c_solver(A(indc, indc), xs); 
                    solver.precondIterations_c = solver.precondIterations_c + iter; 
                else
                    ii = indc; 
                    % agmg(A(ii), x(ii), 20, solver.tolerance, solver.maxIterations, 0, [],- 1); 
                    % agmg(A(ii, ii), x(ii), 20, solver.tolerance, solver.maxIterations, 0, [], 1); 
                    % rs = agmg(A(ii, ii), x(ii), 20, solver.tolerance, solver.maxIterations, 0, [], 3); 
                    % lmax = 20; 
                    % rs = agmg(A(ii, ii), x(ii), lmax*20, 0.01, lmax, 0, [], 0); 
                    rs = A(ii, ii)\xs; 
                end
                ldisp('c end')

            end
            % QUESTION : why is T not included in Gauss-Seidel

            if any(indT)

                ldisp('T start')
                xs = x(indT) - 0.0*A(indT, not(indT))*r(not(indT)); 
                [rT, flag, res, iter] = T_solver(A(indT, indT), xs); 
                solver.precondIterations_T = solver.precondIterations_T + iter; 
                ldisp('T End')
                
                r(indT)   = rT;
            end
            
            r(indphi) = rp; 
            r(indc)   = rs;
            
        end


        function elliptic_solver = getElipticSolver(solver, solverspec) % ok

            switch solverspec.name
                
              case 'agmg'

                nmax = 20;
                elliptic_solver = @(A, b) agmg(A, b, nmax, 1e-5, nmax, 1, [], 0);
                
              case 'agmgsolver_prec'
                % QUESTION : any fundamental difference with 'agmgsolver' above ?
                
                elliptic_solver = @(A, b) solver.agmgprecond(A, b);
                
              case 'direct'
                
                elliptic_solver = @(A, b) deal(mldivide(A,b), 1, 0, 1);
                
              case 'amgcl'
                
                % isolver = struct('type','bicgstab','M',50);

                isolver = struct('type', 'gmres', ...
                                 'M'   , 50);
                
                relaxation=struct('type', 'ilu0');
                %relaxation=struct('type','spai0')
                %relaxation=struct('type',lower(smoother))
                %maxNumCompThreads(np)
                alpha = 0.01;
                
                aggr = struct('eps_strong', alpha);

                coarsening = struct('type'       , 'aggregation', ...
                                    'over_interp', 1.0          , ...
                                    'aggr'       , aggr)

                coarsetarget = 1200;

                coarsening = struct('type'         , 'ruge_stuben', ...
                                    'rs_eps_strong', alpha        , ...
                                    'rs_trunc'     , true         , ...
                                    'rs_eps_trunc' , alpha);
                
                precond = struct('class'        , 'amg'       , ...
                                 'coarsening'   , coarsening  , ...
                                 'relax'        , relaxation  , ...
                                 'coarse_enough', coarsetarget, ...
                                 'max_levels'   , 20          , ...
                                 'ncycle'       , 1           , ...
                                 'npre'         , 1           , ...
                                 'npost'        , 1           , ...
                                 'pre_cycle'    , 0           , ...
                                 'direct_coarse', true);
                
                opts = struct('solver'      , isolver  , ...
                              'precond'     , precond  ,...
                              'reuse_mode'  , 1        , ...
                              'solver_type' , 'regular', ...
                              'write_params', false    , ...
                              'block_size'  , 1        , ...
                              'verbosity'   , 10);
                
                elliptic_solver = @(A, b) solver.amgcl(A, b, opts);
                
              otherwise
                
                error('Wrong solver type for cpr voltage preconditioner')
                
            end
            
        end
        
        function  [x, flag, relres, iter] = amgcl(solver, A, b, opt) % ok

            tol = 1e-5;
            
            [x, extra] =  amgcl(A, b, 'amgcloptions', opt, ...
                                      'blocksize'   , 1  , ...
                                      'tol'         , tol, ...
                                      'maxiter'     , 20);

            flag   = extra.err < tol;
            relres = extra.err;
            iter   = extra.nIter;
            
            disp(['err', num2str(extra.err),' iter ',num2str(extra.nIter)]);
            
        end

        function  [x, flag, relres, iter] = agmgprecond(solver, A, b)% ok
            
            agmg(A, b, 20, solver.tolerance, solver.maxIterations, solver.verbosity, [], -1); 
            agmg(A, b, 20, solver.tolerance, solver.maxIterations, solver.verbosity, [], 1);
            
            [x, flag, relres, iter] = agmg(A, b, 20, 1e-10, 20, 0, [], 3);
            
            iter = 1;
            
        end

        function smoother = getSmoother(smoother_type) % ok
            
            switch smoother_type

              case 'ilu0'
                
                iluopt = struct('type', 'nofill'); % , 'droptol', 0.1); 
                [L, U] = ilu(A(inds, inds), iluopt); 
                opt = struct('L'       , L, ...
                             'U'       , U, ...
                             'smoother', 'ilu0'); 
                dd = abs(diag(U)); 

                if (max(dd)/min(dd) > 1e14)
                    error(); 
                end
                
                smoother = @(A, x) opt.U\(opt.L\x); 

              case 'gs'
               
                Atmp = A(inds, inds); 
                U    = triu(Atmp);
                L    = Atmp - U; 
                dd   = abs(diag(U));
                
                if(max(dd)/min(dd) > 1e14)
                    error(); 
                end
                
                opt = struct('L'       , L   , ...
                             'U'       , U   , ...
                             'ind'     , inds, ...
                             'smoother', 'gs'); 
                rhs = opt.U\x(opt.ind); 
                smoother = @(A, x) - opt.U\(opt.L*(opt.U\x)) + opt.U\x; 
                
              otherwise
                error('smoother_type not recognized');
            end
            
        end
                
        
        function dx = storeIncrements(solver, problem, result) % ok
            % Extract the results from a vector into a cell array with one
            % entry per primary variable in the linearized problem.
            
            % Find first index corresponding to ADI equation
            ix = find(cellfun(@(x) isa(x, 'ADI'), problem.equations), 1);
            if isempty(ix)
                % No ADI equations, we return empty increments
                assert(isempty(result), ...
                    'No ADI equations returned. Unable to map increments to variables.');
                dx = {};
                return
            end
            % Calculate positions in newton increment
            numVars = problem.equations{ix}.getNumVars();
            cumVars = cumsum(numVars);
            ii = [[1;cumVars(1:end-1)+1], cumVars];
            
            eqn = size(ii,1);
            dx = cell(eqn,1);
            for i = 1:eqn
                dx{i} = result(ii(i,1):ii(i,2), :);
            end
        end
        
        function [A, b, sys, x0] = reduceLinearSystem(solver, A, b, isAdjoint, x0) %ok
            % Perform Schur complement reduction of linear system
            if nargin < 4
                isAdjoint = false;
            end
            if nargin < 5
                x0 = [];
            end
            if solver.useSparseReduction
                method = 'sparse';
            else
                method = 'subset';
            end
            sys = splitMatrixForReduction(A, b, solver.keepNumber, method, true);
            % QUESTION : shall we use Splitmatrixforreduction ?
            if ~isempty(sys.B)
                A = sys.B - sys.C*(sys.E_U\(sys.E_L\sys.D));
                if isAdjoint
                    % We are solving the transpose system
                    b = sys.f - (sys.D')*((sys.E')\sys.h);
                else
                    % We are solving the untransposed system
                    b = sys.f - sys.C*(sys.E_U\(sys.E_L\sys.h));
                end
                if ~isempty(x0)
                    x0 = x0(1:solver.keepNumber);
                end
            end
        end

        function x = recoverLinearSystem(solver, x, sys) % ok
            % Recover eliminated variables
            if ~isempty(sys.E_U)
                s = sys.E_U\(sys.E_L\(sys.h - sys.D*x));
                x = [x; s];
            end
        end
        
        function x = recoverLinearSystemAdjoint(solver, x, sys) %ok
            % Recover eliminated variables
            if ~isempty(sys.E)
                s = (sys.E')\(sys.h - sys.C'*x);
                x = [x; s];
            end
        end

        function [A, b, scaling, x0] = applyScaling(solver, A, b, x0) %ok
            % Apply left or right diagonal scaling
            if nargin == 4
                x0 = [];
            end
            scaling = struct();
            applyLeft = solver.applyLeftDiagonalScaling;
            applyRight = solver.applyRightDiagonalScaling;
            if ~applyLeft && ~applyRight
                return
            end
            M = solver.getDiagonalInverse(A);
            if solver.applyLeftDiagonalScaling
                assert(~applyRight, 'Cannot both apply left and right diagonal scaling');
                A = M*A;
                b = M*b;
            else
                A  = A*M;
                if ~isempty(x0)
                    x0 = M\x0;
                end
            end
            scaling.M = M;
        end
        
        function x = undoScaling(solver, x, scaling) %ok
            % Undo effects of scaling applied to linear system
            if solver.applyRightDiagonalScaling
                x = scaling.M*x;
            end
        end
        
        function x = undoScalingAdjoint(solver, x, scaling) % ok
            % Undo effects of scaling applied to linear system (adjoint
            % version)
            if solver.applyLeftDiagonalScaling
                x = scaling.M*x;
            end
        end

        function x = preconditionerInverse(solver, M, x) % keep ?
            % Apply a preconditioner. Either a handle or a matrix.
            if isempty(M)
                return
            end
            
            if isa(M, 'function_handle')
                % We got a function handle for the inverse
                x = M(x);
            else
                % We got a matrix
                x = M\x;
            end
        end
        
        function [A, b, order, x0] = reorderLinearSystem(solver, A, b, order, x0) %ok
        
            if nargin < 5
                x0 = [];
            end
            vo = solver.variableOrdering;
            eo = solver.equationOrdering;
            hasVar = ~isempty(vo);
            hasEq = ~isempty(eo);
            n = size(A, 1);
            if hasVar
                if isa(vo, 'function_handle')
                    vo = vo(A, b);
                end
                nv = numel(vo);
                if nv < n
                    vo = [vo; (nv+1:n)'];
                elseif nv > n
                    vo = vo(1:n);
                end
            end
            if hasEq
                if isa(eo, 'function_handle')
                    eo = eo(A, b);
                elseif isnan(eo)
                    eo = vo;
                end
                ne = numel(eo);
                if ne < n
                    eo = [eo; (ne+1:n)'];
                elseif ne > n
                    eo = eo(1:n);
                end
            end
            if hasVar && hasEq
                A  = A(eo, vo);
                b  = b(eo);
            elseif hasVar
                A  = A(:, vo);
            elseif hasEq
                A = A(eo, :);
                b = b(eo);
            end
            if ~isempty(x0)
                x0 = x0(vo);
            end
            order = struct('variableOrdering', vo, 'equationOrdering', eo);
        end
        
        function x = deorderLinearSystem(solver, x, order) % ok
            if nargin < 3
                tmp = solver;
            else
                tmp = order;
            end
            if ~isempty(solver.variableOrdering)
                x(tmp.variableOrdering) = x(1:numel(tmp.variableOrdering));
            end
        end

        function x = deorderLinearSystemAdjoint(solver, x, order) % ok
            
            if nargin < 3
                tmp = solver;
            else
                tmp = order;
            end
            if ~isempty(solver.equationOrdering)
                x(tmp.equationOrdering) = x(1:numel(tmp.equationOrdering));
            end
            
        end

        function [problem, eliminated] = reduceToVariable(solver, problem, keep) % ok
            
            remove = find(~keep);
            
            problem = problem.clearSystem();
            
            eliminated = cell(numel(remove), 1);
            elimNames = problem.equationNames(remove);
            for i = 1:numel(remove)
                [problem, eliminated{i}] = problem.eliminateVariable(elimNames{i});
            end
            
        end
        
        function dx = recoverResult(solver, dxElim, eliminatedEqs, keep)

            kept = find(keep);
            left = find(~keep);
            nokept = find(~keep);
            keptEqNo = NaN(size(nokept));
            for i=1:numel(nokept)
                keptEqNo(i) = sum(kept<nokept(i));
            end
            
            % Find number of variables
            nP = numel(keep);
            
            % Set up storage for all variables, including those we eliminated previously
            dx = cell(nP, 1);
            
            % Recover non-cell variables
            recovered = false(nP, 1);
            recovered(kept) = true;
            
            % Put the recovered variables into place
            dx(recovered) = dxElim;
            
            for i = numel(eliminatedEqs) : -1 : 1
                pos = left(i);
                dVal = recoverVars(eliminatedEqs{i}, keptEqNo(i) + 1, dx(recovered));
                dx{pos} = dVal;
                recovered(pos) = true;
            end
        end
        
        function M = getDiagonalInverse(solver, A) % ok
            % Reciprocal of diagonal matrix
            sz = size(A);
            assert(sz(1) == sz(2), 'Matrix must be square!');
            n = sz(1);
            d = 1./abs(diag(A));
            d(~isfinite(d)) = 1;
            I = (1:n)';
            M = sparse(I, I, d, n, n);
        end
        
        function initialGuess = getInitialGuess(solver, problem) % ok
            if isempty(solver.initialGuessFun)
                initialGuess = {[]};
                return
            end
            initialGuess = solver.initialGuessFun(problem);
        end
        
        function disp(solver)
            [d, sn] = solver.getDescription();
            s1 = sprintf('  %s linear solver of class %s', sn, class(solver));
            fprintf('%s\n  %s\n  %s\n  ->', s1, repmat('-', 1, numel(s1)-2), d);
            builtin('disp', solver);
        end

        
    end
end





%{
Copyright 2021-2022 SINTEF Industry, Sustainable Energy Technology
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
