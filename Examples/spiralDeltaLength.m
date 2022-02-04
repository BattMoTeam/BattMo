function dlperdx = spiralDeltaLength(x, y, params)

    w  = params.w;
    r0 = params.r0;

    theta = x./r0;
    dlperdx = 1/(2*pi*r0)*sqrt((y + w).^2 + (2*pi*r0 + theta.*(y + w)).^2);
    
end
