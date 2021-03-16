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
                    vars=[2,4,6,7,8,9]
                    % equation 4,6 seems ok.
                    % 2 have very large onditionnumber, 7,8 also seems to
                    % heigh
                    %vars=[4,6]
                    %vars=[]
                    %vars=[7,9];
                    %vars=[8]
                    pos=cumsum(numVars)
                    pos=[[1;pos(1:end-1)+1],pos];
                    posvar=pos(vars,:);
                    ind=mcolon(posvar(:,1),posvar(:,2));
                    indb=false(size(A,1),1);
                    indb(ind) = true;
                    f =@(x) solver.precond(x,A,indb);
                    AA=A(indb,indb);
                    bb=b(indb);
                    %x=agmg(AA,rand(size(bb)),1,[],[],1);
                    %condest(AA)
                    x = gmres(A,b,10,solver.tol,solver.maxiter,f);
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
        
        function r = precond(solver,x,A,ind)
            r=x;
            rphi = A(ind,ind)\x(ind);
            r(ind)=rphi;
        end
    end
end