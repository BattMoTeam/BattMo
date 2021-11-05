classdef LinearSolverBatteryNew < LinearSolverAD
    properties
        method
        verbosity
        first
        reuse_setup
    end
    methods
        function solver = LinearSolverBatteryNew(varargin)            
            opt=struct(...
                'method','direct',...
                'verbosity',0,...
                'reuse_setup',false)
            [opt,extra] = merge_options(opt,varargin{:});
            solver = solver@LinearSolverAD(extra{:});
            solver.method = opt.method;           
            solver.verbosity=opt.verbosity;
            solver.first=true;
            solver.reuse_setup =  opt.reuse_setup;
        end
    
        function [result, report] = solveLinearSystem(solver, A, b, x0, problem)
            report = solver.getSolveReport();
            switch solver.method
                case 'direct'
                    result=A\b;
                case 'agmg'
                    a=tic();
                    if(solver.reuse_setup)                        
                        if(solver.first)
                            solver.first = false;
                            agmg(A,b,20,solver.tolerance,solver.maxIterations,solver.verbosity,[],-1);
                            agmg(A,b,20,solver.tolerance,solver.maxIterations,solver.verbosity,[],1);
                        end
                        [result,flag,relres,iter]=agmg(A,b,20,solver.tolerance,solver.maxIterations,solver.verbosity,[],2);
                        if(flag == 1)
                            agmg(A,b,20,solver.tolerance,solver.maxIterations,solver.verbosity,[],-1);
                            agmg(A,b,20,solver.tolerance,solver.maxIterations,solver.verbosity,[],1);
                            [result,flag,relres,iter_new]=agmg(A,b,20,solver.tolerance,solver.maxIterations,solver.verbosity,result,2);
                            iter = iter+iter_new;
                            solver.first=true;
                        end
                    else
                        [result,flag,relres,iter]=agmg(A,b,20,solver.tolerance,solver.maxIterations,solver.verbosity);
                        
                    end
                    report.Iterations = iter;
                    report.Residual = relres;
                    report.LinearSolutionTime = toc(a);
                    report.Converged = flag; %% should we set always true?                    
                    %if(reset)
                    %    result=agmg(A,b,20,solver.tolerance,solver.maxIterations,solver.verbosity,-1);
                    %end
                case 'matlab_agmg'
                    agmg(A,b,20,solver.tolerance,solver.maxIterations,solver.verbosity,[],-1);
                    agmg(A,b,20,solver.tolerance,solver.maxIterations,solver.verbosity,[],1);
                    f=@(b) solver.amgprecond(b,A)
                    result = gmres(A,b,10,solver.tol,solver.maxiter,f);
                case 'matlab_cpr_amg'
                    indb = getPotentialIndex(problem)
                    agmg(A(indb,indb),b(inb),20,solver.tolerance,solver.maxIterations,solver.verbosity,[],-1);
                    agmg(A(indb,indb),b(indb),20,solver.tolerance,solver.maxIterations,solver.verbosity,[],1);
                    f=@(b) solver.amgprecondcpr(b,A,indb,false)
                    result = gmres(A,b,10,solver.tol,solver.maxiter,f);
                case 'matlab_cpr'
                    indb = solver.getPotentialIndex(problem) ;
                    f=@(b) solver.precondcpr(b,A,indb,true);
                    result = gmres(A,b,20,solver.tolerance,solver.maxIterations,f);
                case 'amgcl'
                      mycase = 'amg';
                      switch mycase
                          case 'ILU0'
                            isolver = struct('type','gmres','M',200,'verbosity',10);
                            precond = struct('class','relaxation','type','ilu0','damping',1,'verbosity',10);
                            options = struct('solver',isolver,'precond',precond,'solver_type','regular',...
                                'write_params',true,'block_size',1,'verbosity',0,'reuse_mode',1);
                          case 'amg'
                              %isolver = struct('type','gmres','M',20,'verbosity',3);
                              isolver = struct('type','bicgstab','M',50);
                               relaxation=struct('type','ilu0')
                                %relaxation=struct('type','spai0')
                                %relaxation=struct('type',lower(smoother))
                                %maxNumCompThreads(np)
                                alpha=0.001;
                                 aggr=struct('eps_strong',alpha);
                                 % coarsening=struct('type','aggregation','over_interp',1.0,'aggr',aggr)
                                 %coarsening=struct('type','smoothed_aggregation','relax',2/3,'aggr',aggr,'estimate_spectral_radius',false,'power_iters',10,'over_interp',1.0)
                                 coarsetarget=1200;
                                 coarsening=struct('type','ruge_stuben','rs_eps_strong',alpha,'rs_trunc',true,'rs_eps_trunc',alpha);
                                 precond = struct('class','amg','coarsening',coarsening,'relax',relaxation,...
                                    'coarse_enough',coarsetarget,'max_levels',20,'ncycle',1,'npre',1,'npost',1,'pre_cycle',0,'direct_coarse',true);
                                 options = struct('solver',isolver,'precond',precond,...
                                     'reuse_mode',1,'solver_type','regular','write_params',false,'block_size',1,'verbosity',10);
                              
                          otherwise
                              error()
                      end
                      tic;
                    [result,extra]=amgcl(A,b,'amgcloptions',options,'blocksize',1,'tol', solver.tolerance,'maxiter',solver.maxIterations);
                    a=struct('iterations',num2str(extra.nIter),'reduction',extra.err)
                    toc;
                otherwise
                    error('Method not implemented');
            end
            
            %% fill report
        end
        
        function indb = getPotentialIndex(solver,problem)
            numVars = problem.equations{1}.getNumVars();
            %eqnames = problem.
            if(numel(problem.equations)>8)
                vars=[2,4,6,7,8,9];
            else
                vars=[2,4,6,7,8];
            end
                    % equation 4,6 seems ok.
            pos=cumsum(numVars);
            pos=[[1;pos(1:end-1)+1],pos];
            posvar=pos(vars,:);
            %eind = mcolon(posvar(1,1),posvar(1,2));
            %neind = mcolon(posvar(2:end,1),posvar(2:end,2));
            ind=mcolon(posvar(:,1),posvar(:,2));
            indb=false(pos(end),1);      
            indb(ind) = true;
        end            
        function r = precondcpr(solver,x,A,ind,use_matlab)%,AA)
            r=x*0;
            %[L,U,p] = lu(A,'vector');
            %% post smooth
            [L,U]=ilu(A);
            dr= U\(L\x);
            %r=A\x;
            r=r+dr;
            x=x-A*dr;
            if(true)
                if(use_matlab)
                    rphi = A(ind,ind)\x(ind);
                else
                    rphi=  agmg(A(ind,ind),x(ind),0,1e-4,20,0,2);
                end
                r(ind)=r(ind)+rphi;
                x(ind)=x(ind)-A(ind,ind)*rphi;
            end
            %% pre smooth
             %% post smooth
            dr= U\(L\x);
            r=r+dr;
            %x=x-A*dr;
        end
        
        function r = agmgprecond(solver,x,A)            
                r= agmg(A,x,0,1e-4,20,0,2);                
        end

        
    end
end