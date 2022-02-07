function [Emid, Tmid, energyDensity, specificEnergy, Energy] = computeEnergyDensity(E, I, T, t, vol, mass)
    dt = diff(t);
    Emid = (E(2:end) + E(1:end - 1))/2.0; 
    Tmid = (T(2:end) + T(1:end - 1))/2.0;
    Imid = (I(2:end) + I(1:end - 1))/2.0; 
    Energy = cumsum((Emid.*Imid).*dt);
    Energy = Energy/hour; % W h
    energyDensity = Energy/vol/1000; %  W h/L
    specificEnergy = Energy/mass; %  W h/kg
end