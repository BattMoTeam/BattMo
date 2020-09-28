function [result] = div(y, Xb)

%size y = N+1, because y is given on cell interfaces
result = diff(y) ./ diff(Xb) ;                      %size: Nx1
