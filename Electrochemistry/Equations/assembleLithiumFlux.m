function state = assembleLithiumFlux(model, state)
    
    D = state.ActiveMaterial.D;
    cLi = state.ActiveMaterial.Li;
    Deff = D .* model.volumeFraction .^1.5;
            
    LiFlux = assembleFlux(model, potential, Deff);
            
    state.LiFlux = flux;
    
end

