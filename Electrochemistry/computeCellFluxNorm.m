function  src = computeCellFluxNorm(model, j, coef)
%
% - Compute cell-valued square of the norm of a face-valued flux (j), weighted by inverse of coefficient (coef)
% - We weight also by volume and volume fraction.
%

    
    op = model.operators.cellFluxOp;
    
    if isprop(model, 'volumeFraction')
        volfrac = model.volumeFraction;
    else
        volfrac = 1;
    end
    vols = model.G.cells.volumes;

    j = op.P*j;
    jsq = j.^2;
    jsq = op.S*jsq;
    
    
    src = (volfrac.*vols).*jsq./coef;
    
    
end