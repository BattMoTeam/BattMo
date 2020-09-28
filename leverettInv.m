function s = leverettInv(par, pc, pf)

    
    J0  = par(1);
    A1  = par(2);
    A2  = par(3);
    B1  = par(4);
    B2  = par(5);

    s = (log(pc ./ pf - J0) - log(A1./A2) + 0.5.*(B1 - B2) ) ./ (B1 - B2);
                      
end