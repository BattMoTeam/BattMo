function D = computeElectrolyteConductivity_Chen2020(c, T)
    
    c = c./1000;
    D = 0.1297 .*c.^3 - 2.51 .*c.^(1.5) + 3.329 .*c;

end



%{
%}
