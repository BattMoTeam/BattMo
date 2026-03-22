function D = computeDiffusionCoefficient_Chen2020(c, T)
    
    c = c./1000;
    D = 8.794 .*10^(-11) .*c .^2 - 3.972 .*10^(-10) .*c + 4.862*10^(-10);

end



%{
%}
