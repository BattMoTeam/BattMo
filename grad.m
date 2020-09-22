function result = grad(y, X, varargin)

[~, n] = size(y);

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
        result = [zeros(1,n); bsxfun(@rdivide, diff(y),  diff(X)); zeros(1,n)]; % gradients on compartment boundaries
        
    elseif strcmpi(conditions,'adiabatic-dirlichet') == 1 
        %% Mixed adiabatic-dirlichet boundry conditions
        % Dirlichet boundary condition at right boundary
        yBCR = varargin{rightBCi + 2};
        % Flux over the right boundary
        dXBCR = X(end) - X(end-1);
        NBCR = (yBCR - y(end)) / dXBCR;

        % Calcualte first derivative
        result = [zeros(1,n); bsxfun(@rdivide, diff(y),  diff(X)); NBCR]; % gradients on compartment boundaries

    elseif strcmpi(conditions,'dirlichet-adiabatic') == 1 
        %% Mixed dirlichet-adiabatic boundry conditions
        % Dirlichet boundary condition at left boundary
        yBCL = varargin{leftBCi + 2};
        % Flux over the right boundary
        dXBCL = X(2) - X(1);
        NBCL = (y(1) - yBCL) / dXBCL;

        % Calcualte first derivative
        result = [NBCL; bsxfun(@rdivide, diff(y),  diff(X)); zeros(1,n)]; % gradients on compartment boundaries

    elseif strcmpi(conditions,'dirlichet-dirlichet') == 1 
        %% Dirlichet-dirlichet boundry conditions
        % Dirlichet boundary condition at left boundary
        yBCL = varargin{leftBCi + 2};
        % Flux over the right boundary
        dXBCL = X(2) - X(1);
        NBCL = (y(1) - yBCL) / dXBCL;

        % Dirlichet boundary condition at right boundary
        yBCR = varargin{rightBCi + 2};
        % Flux over the right boundary
        dXBCR = X(end) - X(end-1);
        NBCR = (yBCR - y(end)) / dXBCR;

        % Calcualte first derivative
        result = [NBCL; bsxfun(@rdivide, diff(y),  diff(X)); NBCR]; % gradients on compartment boundaries

        
end

% % Calcualte first derivative
% result = [  zeros(1,n); ...
%             bsxfun(@rdivide, diff(y),  diff(X)); ...
%             zeros(1,n)]; % gradients on compartment boundaries

end