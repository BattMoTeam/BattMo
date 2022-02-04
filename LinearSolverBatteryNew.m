classdef LinearSolverBatteryNew < LinearSolverAD
    properties
        method
        verbosity
    end
    methods
        function solver = LinearSolverBatteryNew(varargin)            
            opt=struct(...
                'method','direct',...
                'verbosity',0)
            [opt,extra] = merge_options(opt,varargin{:});
            solver = solver@LinearSolverAD(extra{:});
            solver.method = opt.method;           
            solver.verbosity=opt.verbosity;           
        end
    
        function [result, report] = solveLinearSystem(solver, A, b, x0)                           
            switch solver.method    
                case 'direct'
                    result=A\b;
                case 'agmg'
                    result=agmg(A,b,0,solver.tolerance,solver.maxIterations,solver.verbosity);
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
            report = solver.getSolveReport();
            %% fill report
        end
        
        function r = precond(solver,x,A,ind)%,AA)
            r=x*0;
            %[L,U,p] = lu(A,'vector');
            %% post smooth
            [L,U]=ilu(A);
            dr= U\(L\x);
            %r=A\x;
            r=r+dr;
            x=x-A*dr;
            if(true)
                rphi=  agmg(A(ind,ind),x(ind),0,1e-4,20,0);
                r(ind)=r(ind)+rphi;
                x(ind)=x(ind)-A(ind,ind)*rphi;
            end
            %% pre smooth
             %% post smooth
            dr= U\(L\x);
            r=r+dr;
            %x=x-A*dr;
        end
    end
end