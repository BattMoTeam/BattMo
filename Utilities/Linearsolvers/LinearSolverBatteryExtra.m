classdef LinearSolverBatteryExtra < LinearSolverADExtra

    properties
        method
        verbosity
        first
        reuse_setup
        precondIterations_phi
        precondIterations_c
    end
    
    methods
        function solver = LinearSolverBatteryExtra(varargin)
            
            opt=struct('method'     ,'direct',...
                       'verbosity'  ,0       ,...
                       'reuse_setup',false);
            [opt,extra] = merge_options(opt,varargin{:});
            solver = solver@LinearSolverADExtra(extra{:});
            solver.method = opt.method;           
            solver.verbosity=opt.verbosity;
            solver.first=true;
            solver.reuse_setup =  opt.reuse_setup;
            
        end
    
        function [result, report] = solveLinearSystem(solver, A, b, x0, problem)
            
            report = solver.getSolveReport();
            switch solver.method
                case 'direct'                    
                    %indb = getPotentialIndex(solver,problem);
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
                        tic;
                        [result,flag,relres,iter]=agmg(A,b,20,solver.tolerance,solver.maxIterations,solver.verbosity);
                        toc;
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
                    f=@(b) solver.agmgprecond(b,A)
                    a=tic;
                    [result,flags, relres, iter] = gmres(A,b,solver.maxIterations,solver.tolerance,solver.maxIterations,f);
                     report.Iterations = (iter(1)-1)*solver.maxIterations + iter(2);
                    report.Residual = relres;
                    report.Converged = flags;
                    report.LinearSolutionTime = toc(a);
                    %toc;
                case 'matlab_p_gs'
                    solver.precondIterations_phi = 0;
                    solver.precondIterations_c = 0;
                    indb = solver.getPotentialIndex(problem);
                    %agmg(A(indb,indb),b(indb),20,solver.tolerance,solver.maxIterations,solver.verbosity,[],-1);
                    %agmg(A(indb,indb),b(indb),20,solver.tolerance,solver.maxIterations,solver.verbosity,[],1);
                    %e_solver =@(A,b) mldivide(A,b);
                    lverbosity = 0;
                    ltol = 1e-3;
                    lmax = 40;
                    %phi_solver = solver.getElipticSolver('direct');
                    phi_solver = solver.getElipticSolver('agmgsolver');
                    %phi_solver = solver.getElipticSolver('amgclsolver');
                    %c_solver = solver.getElipticSolver('agmgsolver');
                    c_solver = solver.getElipticSolver('amgclsolver');
                    %e_solver =@(A,b) agmg(A,b,lmax,ltol,lmax,lverbosity,[],0);
                    %e_solver =@(A,b) mldivide(A,b);
                    opt = [];
                    f=@(b) solver.p_gs_precond(b,A,indb,phi_solver, c_solver, opt);
                    a=tic;
                    restart = solver.maxIterations;
                    %[result, flags, relres, iter]= gmres(A,b,restart,solver.tolerance,solver.maxIterations,f);
                    %[result, flags, relres, iter]= bicgstab(A,b,solver.tolerance,5,f)
                    [result, flags, relres, iter]= gmres(A,b,restart,solver.tolerance,solver.maxIterations,f);
                    report.Iterations = (iter(1)-1)*solver.maxIterations + iter(2);
                    report.Residual = relres;
                    report.Converged = flags;
                    report.precondIterations_phi = solver.precondIterations_phi; 
                    report.precondIterations_c = solver.precondIterations_c;
                    report.LinearSolutionTime = toc(a);
 
                case 'matlab_cpr_agmg'
                    indb = solver.getPotentialIndex(problem);
                    agmg(A(indb,indb),b(indb),20,solver.tolerance,solver.maxIterations,solver.verbosity,[],-1);
                    agmg(A(indb,indb),b(indb),20,solver.tolerance,solver.maxIterations,solver.verbosity,[],1);
                    %rphi=  agmg(A(indb,indb),b(indb),20,1e-4,20,1,,-1);
                    %rphi=  agmg(A(indb,indb),b(indb),20,1e-4,20,1,[],1);
                    solver_type = 'agmgsolver';
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
                      mycase = 'ILU0';
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
                          case 'amg_cpr'
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
            vars=[];
            for i = 1:numel(problem.primaryVariables)
                pp = problem.primaryVariables{i};
                if strcmp(pp{end},'E') || strcmp(pp{end},'phi')
                    vars =[vars,i];
                end 
            end
            %if(numel(problem.equations)>8)
            %    vars=[2,4,6,7,8,9];
            %else
            %    vars=[2,4,6,7,8];
            %end
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
        function r = precondcpr(solver,x,A,ind,voltage_solver,smoother)%,AA)
            r=x*0;          
            dr = smoother(A,x);
            r=r+dr;
            x=x-A*dr;
            rphi = voltage_solver(A(ind,ind),x(ind));            
            r(ind)=r(ind)+rphi;
            x(ind)=x(ind)-A(ind,ind)*rphi;
            %dr = x;
            dr = smoother(A,x);            
            %% pre smooth
             %% post smooth
            %dr = x;
            %dr(opt.ind)= opt.U\(opt.L\x(opt.ind));
            r=r+dr;
            %x=x-A*dr;
        end
     
         function r = p_gs_precond(solver,x,A,ind,phi_solver,c_solver, opt)%,AA)
            r=x*0;
            %xp agmg(A(ind,ind),x(ind),20,1e-4,20,1,[],2);
            %lverb = false
            xs = zeros(sum(not(ind)),1);
            ldisp =@(tt) disp(tt);
            %ldisp =@(tt) disp('');
            %ldisp('voltage start')
            for kk=1:1
            xp = x(ind) - A(ind,not(ind))*xs;
            if(true)
               [rp,flag, res, iter]  = phi_solver(A(ind,ind),xp);
               solver.precondIterations_phi = solver.precondIterations_phi +iter; 
            else
                ii = ind;
                %agmg(A(ii,ii),xp,20,solver.tolerance,solver.maxIterations,0,[],1);
                rp = agmg(A(ii,ii),xp,20,solver.tolerance,solver.maxIterations,solver.verbosity,[],0);
            end
            ldisp('voltage end')
            xs = x(not(ind)) - 0.0*A(not(ind),ind)*rp;
            ldisp('c start')
            if(true)
                [rs,flag, res, iter] = c_solver(A(not(ind),not(ind)),xs);
                solver.precondIterations_c = solver.precondIterations_c +iter; 
            else
                ii= not(ind);
                %agmg(A(ii),x(ii),20,solver.tolerance,solver.maxIterations,0,[],-1);
                %agmg(A(ii,ii),x(ii),20,solver.tolerance,solver.maxIterations,0,[],1);
                %rs = agmg(A(ii,ii),x(ii),20,solver.tolerance,solver.maxIterations,0,[],3);
                %lmax = 20;
                %rs = agmg(A(ii,ii),x(ii),lmax*20,0.01,lmax,0,[],0);
                rs = A(ii,ii)\xs;
            end
            end
            ldisp('c end')
            r(ind) = rp;
            r(not(ind)) = rs;         
        end
        %function r = agmgprecond(solver,x,A)            
        %        r= agmg(A,x,0,1e-4,20,0,[],3);                
        %end

        function voltage_solver = getElipticSolver(solver,solver_type, opt)
            
            switch solver_type
                case 'agmgsolver'
                    %agmg(A(indb,indb),b(indb),20,solver.tolerance,solver.maxIterations,solver.verbosity,[],-1);
                    %agmg(A(indb,indb),b(indb),20,solver.tolerance,solver.maxIterations,solver.verbosity,[],1);
                    nmax=20;
                    voltage_solver =@(A,b)  agmg(A,b,nmax,1e-5,nmax,0,[],0);
                case 'agmgsolver_prec'
                    %agmg(A,b,20,solver.tolerance,solver.maxIterations,solver.verbosity,[],-1);
                    %agmg(A,b,20,solver.tolerance,solver.maxIterations,solver.verbosity,[],1);
                    voltage_solver =@(A,b)  solver.agmgprecond(A,b);
                case 'direct'
                    voltage_solver =@(A,b) deal(mldivide(A,b),1,0,1);
                case 'amgclsolver'
                    isolver = struct('type','bicgstab','M',50);
                    isolver = struct('type','gmres','M',50);
                    relaxation=struct('type','ilu0');
                    %relaxation=struct('type','spai0')
                    %relaxation=struct('type',lower(smoother))
                    %maxNumCompThreads(np)
                    alpha=0.01;
                    aggr=struct('eps_strong',alpha);
                    coarsening=struct('type','aggregation','over_interp',1.0,'aggr',aggr)
                    %coarsening=struct('type','smoothed_aggregation','relax',2/3,'aggr',aggr,'estimate_spectral_radius',false,'power_iters',10,'over_interp',1.0)
                    coarsetarget=1200;
                    coarsening=struct('type','ruge_stuben','rs_eps_strong',alpha,'rs_trunc',true,'rs_eps_trunc',alpha);
                    precond = struct('class','amg','coarsening',coarsening,'relax',relaxation,...
                        'coarse_enough',coarsetarget,'max_levels',20,'ncycle',1,'npre',1,'npost',1,'pre_cycle',0,'direct_coarse',true);
                    opts = struct('solver',isolver,'precond',precond,...
                        'reuse_mode',1,'solver_type','regular','write_params',false,'block_size',1,'verbosity',10);
                    voltage_solver =@(A,b) solver.amgcl(A,b,opts);
                    %extra
                otherwise
                    error('Wrong solver type for cpr voltage preconditioner')
            end
        end
        
        function  [X,FLAG,RELRES,ITER] = amgcl(solver,A,b,opt)
                    tol = 1e-5
                   [X,extra] =  amgcl(A,b,'amgcloptions', opt,'blocksize', 1,'tol', tol,'maxiter', 20);
                   FLAG = extra.err < tol;
                   RELRES = extra.err;
                   ITER = extra.nIter;
        end

        function  [X,FLAG,RELRES,ITER] = agmgprecond(solver,A,b)
            agmg(A,b,20,solver.tolerance,solver.maxIterations,solver.verbosity,[],-1);
            agmg(A,b,20,solver.tolerance,solver.maxIterations,solver.verbosity,[],1);
            [X,FLAG,RELRES,ITER]  = agmg(A,b,20,1e-10,20,0,[],3);
            ITER = 1;
        end

        function smoother = getSmoother(smoother_type)
            switch smoother_type
                case 'ilu0'
                    iluopt =struct('type','nofill');%,'droptol',0.1);
                    [L,U]=ilu(A(inds,inds),iluopt);
                    opt=struct('L',L,'U',U,'smoother','ilu0');
                    dd=abs(diag(U));
                    if(max(dd)/min(dd)>1e14)
                        error();
                    end
                    smoother =@(A,x) opt.U\(opt.L\x);
                    %gs structure
                case 'gs'
                    Atmp = A(inds,inds);
                    U = triu(Atmp);%- diag(diag(Atmp));
                    L = Atmp-U;
                    dd=abs(diag(U));
                    if(max(dd)/min(dd)>1e14)
                        error();
                    end
                    %invLower = inv(Lower);
                    opt = struct('L',L,'U',U,'ind',inds,'smoother','gs','n');
                    rhs = opt.U\x(opt.ind);
                    %drtmp = 0*rhs;
                    %for i=1:opt.n
                    smoother = @(A,x) -opt.U\(opt.L*(opt.U\x)) + opt.U\x;%%???????????????????????????????????
                    %end
                    %dr(opt.ind) = drtmp;
                otherwise
                    error()
            end
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
