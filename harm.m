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