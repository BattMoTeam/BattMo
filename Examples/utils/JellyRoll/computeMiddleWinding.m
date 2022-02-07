function n = computeMiddleWinding(r0, w, nwindings)

    x  = linspace(0, 2*pi*nwindings*r0, 10000);
    dx = diff(x);

    params.w = w;
    params.r0 = r0;
    dlperdx = spiralDeltaLength(x, 0, params);
    l = cumsum(dlperdx(1 : (end - 1)).*dx);
    
    ind = find(l > l(end)/2, 1, 'first');

    n = ceil(x(ind)/(2*pi*r0));

end