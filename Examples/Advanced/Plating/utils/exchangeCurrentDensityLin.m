function  j0 = exchangeCurrentDensityLin(celyte, celde, T)

    coef = celyte.*celde;
    coef(coef < 0) = 0;
    F = PhysicalConstants.F;
    th = 1;
    j0 = 6.342e-8*0.5*regularizedSqrt(coef, th)*F;
    
end
