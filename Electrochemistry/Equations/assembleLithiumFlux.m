function state = assembleLithiumFlux(model, state)
    
    D = state.ActiveMaterial.D;
    cLi = state.ActiveMaterial.Li;
    
    Deff = D .* model.volumeFraction .^1.5;
            
    trans = op.harmFace(Deff);
    flux = - trans.*op.Grad(cLi);
            
    state.LiFlux = flux;
    
end

