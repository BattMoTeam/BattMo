function y = linreg(x, d)
    
    y = x;
    ind = (x < d);
    if any(ind)
        y(ind) = x(ind).^2;
    end
    
    ind = (x >= d);
    if any(ind)
        y(ind) = x(ind) + d^2 - d;
    end
    
    ind = (x < 0);
    if any(ind)
        y(ind) = 0;
    end
    
end
