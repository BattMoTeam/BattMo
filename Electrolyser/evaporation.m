function r = evaporation(muLiq, muVap, sLiq, kLV, MW)
%EVAPORATION Calculates the rate of evaporation
%   Detailed explanation goes here

if muVap > muLiq
    r = (muLiq - muVap) .* kLV ./ MW;
else
    r = (sLiq) .* (muLiq - muVap) .* kLV ./ MW;
end

end

