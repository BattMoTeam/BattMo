function R = ionomerSorption(kML, aH2OI, aH2OL)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

con =   physicalConstants();

%R   =   (kML ./ MW) .* (con.R .* T .* (log(aH2OI) - log(aH2OL)));
R   =   (kML) .* (aH2OI - aH2OL);
end

