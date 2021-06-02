function y = regularizedSqrt(x, th)

%
%
% SYNOPSIS:
%   function y = regularizedSqrt(x, th)
%
% DESCRIPTION: returns regularized square root by using linear interporation between (0, 0) and (th, sqrt(th))
%
% PARAMETERS:
%   x  - input values
%   th - threshold
%
% RETURNS:
%   y - output values
%

    y = x; % quick way to create y of same dimension as x and also preserved AD
    ind = (x <= th);
    y(~ind) = (x(~ind)).^0.5;
    y(ind) = x(ind)/th*sqrt(th);
    
end

