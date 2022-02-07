function windingnumbers = computeWindingNumbers(fractions, r0, w, nwindings)

    x  = linspace(0, 2*pi*nwindings*r0, 10000);
    dx = diff(x);

    params.w = w;
    params.r0 = r0;
    dlperdx = spiralDeltaLength(x, 0, params);
    l = cumsum(dlperdx(1 : (end - 1)).*dx);
    
    totallength = l(end);
    
    for ind = 1 : numel(fractions)
        xind = find(l/totallength > fractions(ind), 1, 'first');
        windingnumbers(ind) = ceil(x(xind)/(2*pi*r0));
    end
end