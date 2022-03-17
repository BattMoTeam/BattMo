function result = grad(y, X, varargin)
% We have only refactored 'adiabatic-adiabatic' case
[~, n] = size(y);
assert(n == 1, 'this function should only take in vector')
N = numelValue(y); % note that numelValue is overloaded defined for double and AD variables

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

coregrad = 1./diff(X).*(y(2 : N) - y(1 : (N - 1)));
    
if strcmpi(conditions,'adiabatic-adiabatic') == 1
        %% Adiabatic boundry conditions
        % Calculate first derivative
        result = [zeros(1, n); coregrad; zeros(1, n)]; % gradients on compartment boundaries
        
    elseif strcmpi(conditions,'adiabatic-dirichlet') == 1 
        %% Mixed adiabatic-dirichlet boundry conditions
        % Dirichlet boundary condition at right boundary
        yBCR = varargin{rightBCi + 2};
        % Flux over the right boundary
        dXBCR = X(end) - X(end-1);
        NBCR = (yBCR - y(end)) / (dXBCR/2); % ok for uniform mesh

        % Calculate first derivative
        result = [zeros(1,n); coregrad; NBCR]; % gradients on compartment boundaries

    elseif strcmpi(conditions,'dirichlet-adiabatic') == 1 
        %% Mixed dirichlet-adiabatic boundry conditions
        % Dirichlet boundary condition at left boundary
        yBCL = varargin{leftBCi + 2};
        % Flux over the right boundary
        dXBCL = X(2) - X(1);
        NBCL = (y(1) - yBCL) / (dXBCL/2); % ok for uniform mesh

        % Calculate first derivative
        result = [NBCL; coregrad; zeros(1,n)]; % gradients on compartment boundaries

    elseif strcmpi(conditions,'dirichlet-dirichlet') == 1 
        %% Dirichlet-dirichlet boundry conditions
        % Dirichlet boundary condition at left boundary
        yBCL = varargin{leftBCi + 2};
        % Flux over the right boundary
        dXBCL = X(2) - X(1);
        NBCL = (y(1) - yBCL) / (dXBCL/2); % ok for uniform mesh

        % Dirichlet boundary condition at right boundary
        yBCR = varargin{rightBCi + 2};
        % Flux over the right boundary
        dXBCR = X(end) - X(end-1);
        NBCR = (yBCR - y(end)) / (dXBCR/2); % ok for uniform mesh

        % Calculate first derivative
        result = [NBCL; coregrad; NBCR]; % gradients on compartment boundaries

        
end

% % Calculate first derivative
% result = [  zeros(1,n); ...
%             bsxfun(@rdivide, diff(y),  diff(X)); ...
%             zeros(1,n)]; % gradients on compartment boundaries

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
