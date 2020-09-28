function r = evaporation(muLiq, muVap, sLiq, kLV, MW)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if muVap > muLiq
    r = (muLiq - muVap) .* kLV ./ MW;
else
    r = (sLiq) .* (muLiq - muVap) .* kLV ./ MW;
end

end

