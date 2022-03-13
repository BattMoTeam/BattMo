function res = harm(A, X, Xb)

N = length(X);

[~,n] = size(A);

B = (Xb(2:N)- X(1:N-1))./diff(X);

% res = [zeros(1,n); ...
%     A(1:N-1,:).*A(2:N,:)./(bsxfun(@times,B,A(1:N-1,:)) + bsxfun(@times,(1-B),A(2:N,:))); ...
%     zeros(1,n)]; % harmonic mean of A on compartment boundaries

res1 = B./A(1 : (N - 1));
res2 = (1 - B)./A(2 : N);
res = [A(1); 1./(res1 + res2); A(N)]; % harmonic mean of A on compartment boundaries

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
