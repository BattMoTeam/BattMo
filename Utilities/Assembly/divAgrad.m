function res = divAgrad(y, A, X, Xb, varargin)
%DIVAGRAD Calculates the divergence of the product of a
%gradient. 
%   res = divAgrad(obj, y, A, domain) calculates the divergence
%   of the product of the constant A and the gradient of y, res
%   = div(A*dy/dx), in the specified domain.  

N = length(X);
conditions = [];

leftBCi = find(strcmpi(varargin, 'left'), 1);
rightBCi = find(strcmpi(varargin, 'right'), 1);

    if isempty(varargin) == 1
        conditions = 'adiabatic-adiabatic';
    elseif isempty(leftBCi) == 1 && isempty(rightBCi) == 0
        if strcmpi(varargin(rightBCi+1), 'dirichlet') == 1
            conditions = 'adiabatic-dirichlet';
        elseif strcmpi(varargin(rightBCi+1), 'neumann') == 1
            conditions = 'adiabatic-neumann';
        else
            warning('Unrecognized boundary condition. Please define boundary conditions as either dirichlet or neumann.')
        end
    elseif isempty(leftBCi) == 0 && isempty(rightBCi) == 1
        if strcmpi(varargin(leftBCi+1), 'dirichlet') == 1
            conditions = 'dirichlet-adiabatic';
        elseif strcmpi(varargin(leftBCi+1), 'neumann') == 1
            conditions = 'neumann-adiabatic';
        else
            warning('Unrecognized boundary condition. Please define boundary conditions as either dirichlet or neumann.')
        end
    elseif isempty(leftBCi) == 0 && isempty(rightBCi) == 0
        if strcmpi(varargin(leftBCi+1), 'dirichlet') == 1 && strcmpi(varargin(rightBCi+1), 'dirichlet') == 1
            conditions = 'dirichlet-dirichlet';
        elseif strcmpi(varargin(leftBCi+1), 'dirichlet') == 1 && strcmpi(varargin(rightBCi+1), 'neumann') == 1
            conditions = 'dirichlet-neumann';
        elseif strcmpi(varargin(leftBCi+1), 'neumann') == 1 && strcmpi(varargin(rightBCi+1), 'dirichlet') == 1
            conditions = 'neumann-dirichlet';
        elseif strcmpi(varargin(leftBCi+1), 'neumann') == 1 && strcmpi(varargin(rightBCi+1), 'neumann') == 1
            conditions = 'neumann-neumann';
        else
            warning('Unrecognized boundary condition. Please define boundary conditions as either dirichlet or neumann.')
        end
    end

    coregrad = 1./diff(X).*loc_diff(y);
    
    if strcmpi(conditions,'adiabatic-adiabatic') == 1
        %% Adiabatic boundry conditions
        % Calcualte first derivative
        grad_boundary = [0; coregrad; 0]; % gradients on compartment boundaries
        B = (Xb(2:N)- X(1:N-1))./diff(X);
        A_boundary = [0; A(1:N-1).*A(2:N)./(B.*A(1:N-1) + (1-B).*A(2:N)); 0]; % harmonic mean of A on compartment boundaries
        
        % Calcualte second derivative
        Agrad = A_boundary.*grad_boundary;
        
        res = loc_diff(Agrad) ./ diff(Xb); % second derivative
        
    elseif strcmpi(conditions,'adiabatic-dirichlet') == 1 
        %% Mixed adiabatic-dirichlet boundry conditions
        % Dirichlet boundary condition at right boundary
        yBCR = varargin{rightBCi + 2};
        % Flux over the right boundary
        dXBCR = X(end) - X(end-1);
        NBCR = (yBCR - y(end)) / dXBCR;
        % Value of A at right boundary
        ABCR = A(end);%interp1(X,A,X(end)+0.5*dXBCR,'spline','extrap');

        % Calcualte first derivative
        grad_boundary = [0; coregrad; NBCR]; % gradients on compartment boundaries
        B = (Xb(2:N)- X(1:N-1))./diff(X);
        A_boundary = [0; A(1:N-1).*A(2:N)./(B.*A(1:N-1) + (1-B).*A(2:N)); ABCR]; % harmonic mean of A on compartment boundaries

        % Calcualte second derivative
        Agrad = A_boundary.*grad_boundary;
        res = loc_diff(Agrad) ./ diff(Xb); % second derivative

    elseif strcmpi(conditions,'adiabatic-neumann') == 1 
        %% Mixed adiabatic-neumann boundry conditions
        % Flux over the right boundary
        NBCR = varargin{rightBCi + 2};

        % Calcualte first derivative
        grad_boundary = [0; coregrad; NBCR]; % gradients on compartment boundaries
        B = (Xb(2:N)- X(1:N-1))./diff(X);
        A_boundary = [0; A(1:N-1).*A(2:N)./(B.*A(1:N-1) + (1-B).*A(2:N)); ABCR]; % harmonic mean of A on compartment boundaries

        % Calcualte second derivative
        Agrad = A_boundary.*grad_boundary;
        Agrad(end) = NBCR;
        res = loc_diff(Agrad) ./ diff(Xb); % second derivative

    elseif strcmpi(conditions,'dirichlet-adiabatic') == 1 
        %% Mixed dirichlet-adiabatic boundry conditions
        % Dirichlet boundary condition at left boundary
        yBCL = varargin{leftBCi + 2};
        % Flux over the right boundary
        dXBCL = X(2) - X(1);
        NBCL = (y(1) - yBCL) / dXBCL;
        % Value of A at right boundary
        ABCL = A(1); %interp1(X,A,X(1)-0.5*dXBCL,'spline','extrap');

        % Calcualte first derivative
        grad_boundary = [NBCL; coregrad; 0]; % gradients on compartment boundaries
        B = (Xb(2:N)- X(1:N-1))./diff(X);
        A_boundary = [ABCL; A(1:N-1).*A(2:N)./(B.*A(1:N-1) + (1-B).*A(2:N)); 0]; % harmonic mean of A on compartment boundaries

        % Calcualte second derivative
        Agrad = A_boundary.*grad_boundary;
        res = loc_diff(Agrad) ./ diff(Xb); % second derivative

    elseif strcmpi(conditions,'neumann-adiabatic') == 1 
        %% Mixed neumann-adiabatic boundry conditions
        % Flux over the left boundary
        NBCL = varargin{leftBCi + 2};

        % Calcualte first derivative
        grad_boundary = [0; coregrad; 0]; % gradients on compartment boundaries
        B = (Xb(2:N)- X(1:N-1))./diff(X);
        A_boundary = [0; A(1:N-1).*A(2:N)./(B.*A(1:N-1) + (1-B).*A(2:N)); 0]; % harmonic mean of A on compartment boundaries

        % Calcualte second derivative
        Agrad = A_boundary.*grad_boundary;
        Agrad(1) = NBCL;
        res = loc_diff(Agrad) ./ diff(Xb); % second derivative

    elseif strcmpi(conditions,'dirichlet-neumann') == 1 
        %% Mixed dirichlet-neumann boundry conditions
        % Dirichlet boundary condition at left boundary
        yBCL = varargin{leftBCi + 2};
        % Flux over the right boundary
        dXBCL = X(2) - X(1);
        NBCL = (y(1) - yBCL) / dXBCL;
        % Value of A at right boundary
        ABCL = interp1(X,A,X(1)-0.5*dXBCL,'spline','extrap');

        % Flux over the right boundary
        NBCR = varargin{rightBCi + 2};

        % Calcualte first derivative
        grad_boundary = [NBCL; coregrad; 0]; % gradients on compartment boundaries
        B = (Xb(2:N)- X(1:N-1))./diff(X);
        A_boundary = [ABCL; A(1:N-1).*A(2:N)./(B.*A(1:N-1) + (1-B).*A(2:N)); 0]; % harmonic mean of A on compartment boundaries

        % Calcualte second derivative
        Agrad = A_boundary.*grad_boundary;
        Agrad(end) = NBCR;
        res = loc_diff(Agrad) ./ diff(Xb); % second derivative

    elseif strcmpi(conditions,'dirichlet-dirichlet') == 1 
        %% Dirichlet-dirichlet boundry conditions
        % Dirichlet boundary condition at left boundary
        yBCL = varargin{leftBCi + 2};
        % Flux over the right boundary
        dXBCL = X(2) - X(1);
        NBCL = (y(1) - yBCL) / dXBCL;
        % Value of A at right boundary
        ABCL = interp1(X,A,X(1)-0.5*dXBCL,'spline','extrap');

        % Dirichlet boundary condition at right boundary
        yBCR = varargin{rightBCi + 2};
        % Flux over the right boundary
        dXBCR = X(end) - X(end-1);
        NBCR = (yBCR - y(end)) / dXBCR;
        % Value of A at right boundary
        ABCR = interp1(X,A,X(end)+0.5*dXBCR,'spline','extrap');

        % Calcualte first derivative
        grad_boundary = [NBCL; coregrad; NBCR]; % gradients on compartment boundaries
        B = (Xb(2:N)- X(1:N-1))./diff(X);
        A_boundary = [ABCL; A(1:N-1).*A(2:N)./(B.*A(1:N-1) + (1-B).*A(2:N)); ABCR]; % harmonic mean of A on compartment boundaries

        % Calcualte second derivative
        Agrad = A_boundary.*grad_boundary;
        res = loc_diff(Agrad) ./ diff(Xb); % second derivative

    elseif strcmpi(conditions,'neumann-dirichlet') == 1 
        %% Mixed neumann-dirichlet boundry conditions
        % Flux over the left boundary
        NBCL = varargin{leftBCi + 2};

        % Dirichlet boundary condition at right boundary
        yBCR = varargin{rightBCi + 2};
        % Flux over the right boundary
        dXBCR = X(end) - X(end-1);
        NBCR = (yBCR - y(end)) / dXBCR;
        % Value of A at right boundary
        ABCR = interp1(X,A,X(end)+0.5*dXBCR,'spline','extrap');

        % Calcualte first derivative
        grad_boundary = [0; coregrad; NBCR]; % gradients on compartment boundaries
        B = (Xb(2:N)- X(1:N-1))./diff(X);
        A_boundary = [0; A(1:N-1).*A(2:N)./(B.*A(1:N-1) + (1-B).*A(2:N)); ABCR]; % harmonic mean of A on compartment boundaries

        % Calcualte second derivative
        Agrad = A_boundary.*grad_boundary;
        Agrad(1) = NBCL;
        res = loc_diff(Agrad) ./ diff(Xb); % second derivative

    elseif strcmpi(conditions,'neumann-neumann') == 1 
        %% Neumann-neumann boundry conditions
        % Flux over the left boundary
        NBCL = varargin{leftBCi + 2};

        % Flux over the right boundary
        NBCR = varargin{rightBCi + 2};

        % Calcualte first derivative
        grad_boundary = [0; coregrad; 0]; % gradients on compartment boundaries
        B = (Xb(2:N)- X(1:N-1))./diff(X);
        A_boundary = [0; A(1:N-1).*A(2:N)./(B.*A(1:N-1) + (1-B).*A(2:N)); 0]; % harmonic mean of A on compartment boundaries

        % Calcualte second derivative
        Agrad = A_boundary.*grad_boundary;
        Agrad(1) = NBCL;
        Agrad(end) = NBCR;
        res = loc_diff(Agrad) ./ diff(Xb); % second derivative

    end

end

function dy = loc_diff(y)
% variant of diff that is compatible for AD
    N = numelValue(y);
    dy = y(2 : N) - y(1 : (N - 1));
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
