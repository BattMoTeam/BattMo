function conductivity = computeElectrolyteConductivity_Nyman2008(c, T)
    
    conductivity = 0.1297*(c./1000).^3 - 2.51*(c./1000).^1.5 + 3.329*(c./1000);
    
end


%{
%}
