function state = assembleLithiumFlux(model, state)
    
    D = state.ActiveMaterial.D;
    c = state.ActiveMaterial.Li;
    Deff = D .* model.volumeFraction .^1.5;
            
    LiFlux = assembleFlux(model, c, Deff);
            
    state.LiFlux = LiFlux;
    
end

