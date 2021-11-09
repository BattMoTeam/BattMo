function  src = computeCellFluxNorm(model, j, volumeFraction, coefficient)
%
% - Compute cell-valued square of the norm of a face-valued flux (j), weighted by inverse of coefficient (coef)
% - We weight also by volume and volume fraction.
%

    
    op = model.operators.cellFluxOp;

    vols = model.G.cells.volumes;

    j = op.P*j;
    jsq = j.^2;
    jsq = op.S*jsq;
    
    src = (volumeFraction.*vols).*jsq./coefficient;
    
    
end