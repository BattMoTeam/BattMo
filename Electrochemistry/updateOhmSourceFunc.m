function  state = updateOhmSourceFunc(model, state)

    op = model.operators;
    R = model.ohmicResistance;
    vols = G.cells.volumes;
    
    j = state.j;
    jsq = op.faceToNorm(j);
    jsq = jsq.^2;
    
    state.jHeatOhmSource = vols.*(jsq.^2)./R;
    
end