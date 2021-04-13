classdef LinearSolverBattery 
    properties
        method
        maxiter
        tol
    end
    methods
        function solver = LinearSolverBattery(varargin)
            opt=struct(...
                'method','direct',...
                'tol',1e-4,'maxiter',40);
            opt = merge_options(opt,varargin{:});
            solver.method = opt.method;
            solver.maxiter = opt.maxiter;
            solver.tol = opt.tol;
        end
    
        function [dx, result, report] = solveLinearProblem(solver, problem, model)
            [A, b] = problem.getLinearSystem();                        
            switch solver.method    
                case 'direct'
                    result=A\b;
                case 'iterative'
                    %%
                    numVars = problem.equations{1}.getNumVars();%
                    vars=[2,4,6,7,8,9];
                    % equation 4,6 seems ok.
                    % 2 have very large onditionnumber, 7,8 also seems to
                    % heigh
                    %vars=[4,6]
                    %vars=[]
                    %vars=[7,9];
                    %vars=[8]
                    pos=cumsum(numVars);
                    pos=[[1;pos(1:end-1)+1],pos];
                    posvar=pos(vars,:);
                    eind = mcolon(posvar(1,1),posvar(1,2));
                    neind = mcolon(posvar(2:end,1),posvar(2:end,2));
                    ind=mcolon(posvar(:,1),posvar(:,2));
                    indb=false(size(A,1),1);                    
                    indb(ind) = true;
                    
                    AA=A(indb,indb);
                    bb=b(indb);
                    dv=posvar(:,2)-posvar(:,1);pos=cumsum(dv+1);ptmp=[[0;pos(1:end-1)]+1,pos(1:end)];
                    %eind = mcolon(ptmp(1,1),ptmp(1,2));
                    %neind = mcolon(ptmp(2:end,1),ptmp(2:end,2));
                    %AAtmp=0*AA;
                    %AAtmp(eind,neind)=AA(eind,neind);
                    %AAtmp(neind,eind)=AA(neind,eind);
                    %AAtmp = AAtmp - diag(sum(AAtmp,2));
                    %AAorg=AA;
                    %AA=AA-AAtmp;
                    f =@(x) solver.precond(x,A,indb);%,AA);
                    %x=agmg(AA,rand(size(bb)),1,[],[],1);
                    %condest(AA)
                    result = gmres(A,b,10,solver.tol,solver.maxiter,f);
                    %eig(full(AA));
                    %max(abs(x-result))
                    %A(1,:)
                    %%
                otherwise
                    error('Method not implemented');
            end
            report = [];
            numVars = problem.equations{1}.getNumVars();% depend on and
            cumVars = cumsum(numVars);
            ii = [[1;cumVars(1:end-1)+1], cumVars];
            
            eqn = size(ii,1);
            dx = cell(eqn,1);
            for i = 1:eqn
                dx{i} = result(ii(i,1):ii(i,2), :);
            end
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
            %rphi = A(ind,ind)\x(ind);
            %rphi = AA\x(ind);
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