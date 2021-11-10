function  jsq = computeCellFluxNorm(model, j)
%
% - Compute cell-valued square of the norm of a face-valued flux (j)
%
    
    op = model.operators.cellFluxOp;

    j = op.P*j;
    jsq = j.^2;
    jsq = op.S*jsq;
    
end