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
        if strcmpi(varargin(rightBCi+1), 'dirlichet') == 1
            conditions = 'adiabatic-dirlichet';
        elseif strcmpi(varargin(rightBCi+1), 'neumann') == 1
            conditions = 'adiabatic-neumann';
        else
            warning('Unrecognized boundary condition. Please define boundary conditions as either dirlichet or neumann.')
        end
    elseif isempty(leftBCi) == 0 && isempty(rightBCi) == 1
        if strcmpi(varargin(leftBCi+1), 'dirlichet') == 1
            conditions = 'dirlichet-adiabatic';
        elseif strcmpi(varargin(leftBCi+1), 'neumann') == 1
            conditions = 'neumann-adiabatic';
        else
            warning('Unrecognized boundary condition. Please define boundary conditions as either dirlichet or neumann.')
        end
    elseif isempty(leftBCi) == 0 && isempty(rightBCi) == 0
        if strcmpi(varargin(leftBCi+1), 'dirlichet') == 1 && strcmpi(varargin(rightBCi+1), 'dirlichet') == 1
            conditions = 'dirlichet-dirlichet';
        elseif strcmpi(varargin(leftBCi+1), 'dirlichet') == 1 && strcmpi(varargin(rightBCi+1), 'neumann') == 1
            conditions = 'dirlichet-neumann';
        elseif strcmpi(varargin(leftBCi+1), 'neumann') == 1 && strcmpi(varargin(rightBCi+1), 'dirlichet') == 1
            conditions = 'neumann-dirlichet';
        elseif strcmpi(varargin(leftBCi+1), 'neumann') == 1 && strcmpi(varargin(rightBCi+1), 'neumann') == 1
            conditions = 'neumann-neumann';
        else
            warning('Unrecognized boundary condition. Please define boundary conditions as either dirlichet or neumann.')
        end
    end

    if strcmpi(conditions,'adiabatic-adiabatic') == 1
        %% Adiabatic boundry conditions
        % Calcualte first derivative
        grad_boundary = [0; diff(y) ./ diff(X); 0]; % gradients on compartment boundaries
        B = (Xb(2:N)- X(1:N-1))./diff(X);
        A_boundary = [0; A(1:N-1).*A(2:N)./(B.*A(1:N-1) + (1-B).*A(2:N)); 0]; % harmonic mean of A on compartment boundaries
        
        % Calcualte second derivative
        Agrad = A_boundary.*grad_boundary;
        res = diff(Agrad) ./ diff(Xb); % second derivative
    elseif strcmpi(conditions,'adiabatic-dirlichet') == 1 
        %% Mixed adiabatic-dirlichet boundry conditions
        % Dirlichet boundary condition at right boundary
        yBCR = varargin{rightBCi + 2};
        % Flux over the right boundary
        dXBCR = X(end) - X(end-1);
        NBCR = (yBCR - y(end)) / dXBCR;
        % Value of A at right boundary
        ABCR = A(end);%interp1(X,A,X(end)+0.5*dXBCR,'spline','extrap');

        % Calcualte first derivative
        grad_boundary = [0; diff(y) ./ diff(X); NBCR]; % gradients on compartment boundaries
        B = (Xb(2:N)- X(1:N-1))./diff(X);
        A_boundary = [0; A(1:N-1).*A(2:N)./(B.*A(1:N-1) + (1-B).*A(2:N)); ABCR]; % harmonic mean of A on compartment boundaries

        % Calcualte second derivative
        Agrad = A_boundary.*grad_boundary;
        res = diff(Agrad) ./ diff(Xb); % second derivative

    elseif strcmpi(conditions,'adiabatic-neumann') == 1 
        %% Mixed adiabatic-neumann boundry conditions
        % Flux over the right boundary
        NBCR = varargin{rightBCi + 2};

        % Calcualte first derivative
        grad_boundary = [0; diff(y) ./ diff(X); NBCR]; % gradients on compartment boundaries
        B = (Xb(2:N)- X(1:N-1))./diff(X);
        A_boundary = [0; A(1:N-1).*A(2:N)./(B.*A(1:N-1) + (1-B).*A(2:N)); ABCR]; % harmonic mean of A on compartment boundaries

        % Calcualte second derivative
        Agrad = A_boundary.*grad_boundary;
        Agrad(end) = NBCR;
        res = diff(Agrad) ./ diff(Xb); % second derivative

    elseif strcmpi(conditions,'dirlichet-adiabatic') == 1 
        %% Mixed dirlichet-adiabatic boundry conditions
        % Dirlichet boundary condition at left boundary
        yBCL = varargin{leftBCi + 2};
        % Flux over the right boundary
        dXBCL = X(2) - X(1);
        NBCL = (y(1) - yBCL) / dXBCL;
        % Value of A at right boundary
        ABCL = A(1); %interp1(X,A,X(1)-0.5*dXBCL,'spline','extrap');

        % Calcualte first derivative
        grad_boundary = [NBCL; diff(y) ./ diff(X); 0]; % gradients on compartment boundaries
        B = (Xb(2:N)- X(1:N-1))./diff(X);
        A_boundary = [ABCL; A(1:N-1).*A(2:N)./(B.*A(1:N-1) + (1-B).*A(2:N)); 0]; % harmonic mean of A on compartment boundaries

        % Calcualte second derivative
        Agrad = A_boundary.*grad_boundary;
        res = diff(Agrad) ./ diff(Xb); % second derivative

    elseif strcmpi(conditions,'neumann-adiabatic') == 1 
        %% Mixed neumann-adiabatic boundry conditions
        % Flux over the left boundary
        NBCL = varargin{leftBCi + 2};

        % Calcualte first derivative
        grad_boundary = [0; diff(y) ./ diff(X); 0]; % gradients on compartment boundaries
        B = (Xb(2:N)- X(1:N-1))./diff(X);
        A_boundary = [0; A(1:N-1).*A(2:N)./(B.*A(1:N-1) + (1-B).*A(2:N)); 0]; % harmonic mean of A on compartment boundaries

        % Calcualte second derivative
        Agrad = A_boundary.*grad_boundary;
        Agrad(1) = NBCL;
        res = diff(Agrad) ./ diff(Xb); % second derivative

    elseif strcmpi(conditions,'dirlichet-neumann') == 1 
        %% Mixed dirlichet-neumann boundry conditions
        % Dirlichet boundary condition at left boundary
        yBCL = varargin{leftBCi + 2};
        % Flux over the right boundary
        dXBCL = X(2) - X(1);
        NBCL = (y(1) - yBCL) / dXBCL;
        % Value of A at right boundary
        ABCL = interp1(X,A,X(1)-0.5*dXBCL,'spline','extrap');

        % Flux over the right boundary
        NBCR = varargin{rightBCi + 2};

        % Calcualte first derivative
        grad_boundary = [NBCL; diff(y) ./ diff(X); 0]; % gradients on compartment boundaries
        B = (Xb(2:N)- X(1:N-1))./diff(X);
        A_boundary = [ABCL; A(1:N-1).*A(2:N)./(B.*A(1:N-1) + (1-B).*A(2:N)); 0]; % harmonic mean of A on compartment boundaries

        % Calcualte second derivative
        Agrad = A_boundary.*grad_boundary;
        Agrad(end) = NBCR;
        res = diff(Agrad) ./ diff(Xb); % second derivative

    elseif strcmpi(conditions,'dirlichet-dirlichet') == 1 
        %% Dirlichet-dirlichet boundry conditions
        % Dirlichet boundary condition at left boundary
        yBCL = varargin{leftBCi + 2};
        % Flux over the right boundary
        dXBCL = X(2) - X(1);
        NBCL = (y(1) - yBCL) / dXBCL;
        % Value of A at right boundary
        ABCL = interp1(X,A,X(1)-0.5*dXBCL,'spline','extrap');

        % Dirlichet boundary condition at right boundary
        yBCR = varargin{rightBCi + 2};
        % Flux over the right boundary
        dXBCR = X(end) - X(end-1);
        NBCR = (yBCR - y(end)) / dXBCR;
        % Value of A at right boundary
        ABCR = interp1(X,A,X(end)+0.5*dXBCR,'spline','extrap');

        % Calcualte first derivative
        grad_boundary = [NBCL; diff(y) ./ diff(X); NBCR]; % gradients on compartment boundaries
        B = (Xb(2:N)- X(1:N-1))./diff(X);
        A_boundary = [ABCL; A(1:N-1).*A(2:N)./(B.*A(1:N-1) + (1-B).*A(2:N)); ABCR]; % harmonic mean of A on compartment boundaries

        % Calcualte second derivative
        Agrad = A_boundary.*grad_boundary;
        res = diff(Agrad) ./ diff(Xb); % second derivative

    elseif strcmpi(conditions,'neumann-dirlichet') == 1 
        %% Mixed neumann-dirlichet boundry conditions
        % Flux over the left boundary
        NBCL = varargin{leftBCi + 2};

        % Dirlichet boundary condition at right boundary
        yBCR = varargin{rightBCi + 2};
        % Flux over the right boundary
        dXBCR = X(end) - X(end-1);
        NBCR = (yBCR - y(end)) / dXBCR;
        % Value of A at right boundary
        ABCR = interp1(X,A,X(end)+0.5*dXBCR,'spline','extrap');

        % Calcualte first derivative
        grad_boundary = [0; diff(y) ./ diff(X); NBCR]; % gradients on compartment boundaries
        B = (Xb(2:N)- X(1:N-1))./diff(X);
        A_boundary = [0; A(1:N-1).*A(2:N)./(B.*A(1:N-1) + (1-B).*A(2:N)); ABCR]; % harmonic mean of A on compartment boundaries

        % Calcualte second derivative
        Agrad = A_boundary.*grad_boundary;
        Agrad(1) = NBCL;
        res = diff(Agrad) ./ diff(Xb); % second derivative

    elseif strcmpi(conditions,'neumann-neumann') == 1 
        %% Neumann-neumann boundry conditions
        % Flux over the left boundary
        NBCL = varargin{leftBCi + 2};

        % Flux over the right boundary
        NBCR = varargin{rightBCi + 2};

        % Calcualte first derivative
        grad_boundary = [0; diff(y) ./ diff(X); 0]; % gradients on compartment boundaries
        B = (Xb(2:N)- X(1:N-1))./diff(X);
        A_boundary = [0; A(1:N-1).*A(2:N)./(B.*A(1:N-1) + (1-B).*A(2:N)); 0]; % harmonic mean of A on compartment boundaries

        % Calcualte second derivative
        Agrad = A_boundary.*grad_boundary;
        Agrad(1) = NBCL;
        Agrad(end) = NBCR;
        res = diff(Agrad) ./ diff(Xb); % second derivative

    end

end

