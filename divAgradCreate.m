function [res,Trans] = divAgradCreate(y, A, X, Xb, varargin)
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

    coregrad = 1./diff(X).*loc_diff(y);
    
    if strcmpi(conditions,'adiabatic-adiabatic') == 1
        %% Adiabatic boundry conditions
        % Calcualte first derivative
        grad_boundary = [0; coregrad; 0]; % gradients on compartment boundaries
        B = (Xb(2:N)- X(1:N-1))./diff(X);
        A_boundary = [0; A(1:N-1).*A(2:N)./(B.*A(1:N-1) + (1-B).*A(2:N)); 0]; % harmonic mean of A on compartment boundaries
        
        Trans = A(1:N-1).*A(2:N)./(B.*A(1:N-1) + (1-B).*A(2:N))./diff(X);
        % Calcualte second derivative
        Agrad = A_boundary.*grad_boundary;
        
        res = loc_diff(Agrad) ./ diff(Xb); % second derivative
        
    else
        error();
    end

end

function dy = loc_diff(y)
% variant of diff that is compatible for AD
    N = numelValue(y);
    dy = y(2 : N) - y(1 : (N - 1));
end
