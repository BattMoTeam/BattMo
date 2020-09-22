function res = sineup(y1, y2, x1, x2, x)
    %SINEUP Creates a sine ramp function
    %
    %   res = sineup(y1, y2, x1, x2, x) creates a sine ramp
    %   starting at value y1 at point x1 and ramping to value y2 at
    %   point x2 over the vector x.
    
    dy = y1-y2;
    dx = abs(x1-x2);

    res = (x >= x1 & x<= x2) .* (dy/2.*cos(pi.*(x - x1)./dx) + y1-(dy/2)) + ...
        (x>x2) .* y2 + ...
        (x<x1) .* y1;
    
end