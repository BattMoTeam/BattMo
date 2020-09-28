function [result] = div(y, Xb)

% size y = N+1, because y is given on cell interfaces
    N = numval(y);
    result = (y(2 : N) - y(1 : (N - 1))) ./ diff(Xb) ;                      %size: Nx1

end
