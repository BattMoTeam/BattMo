function R = ionomerSorption(kML, aH2OI, aH2OL)
%%  TODO : add comment
    warning('check this function (looks weird)!')
    con =   physicalConstants();

    % R   =   (kML ./ MW) .* (con.R .* T .* (log(aH2OI) - log(aH2OL)));
    R   =   (kML) .* (aH2OI - aH2OL);
    
end

