function dlperdx = spiralDeltaLength(x, y, params)

    w  = params.w;
    rInner = params.rInner;

    theta = x./rInner;
    dlperdx = 1/(2*pi*rInner)*sqrt((y + w).^2 + (2*pi*rInner + theta.*(y + w)).^2);
    
end
