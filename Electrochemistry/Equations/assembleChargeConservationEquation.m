function chargeCons = assembleChargeConservationEquation(model, phi, jBcSoure, eSource)
    
    op = model.operators;
    sigmaeff = model.EffectiveElectronicConductivity;
    
    chargeCont = (op.Div(j) - jBcSource)./ model.G.cells.volumes./model.constants.F - eSource;
    
    j = - op.harmFace(sigmaeff).*op.Grad(phi); 
    
    chargeCons = (op.Div(j) - jBcSource)./ model.G.cells.volumes./model.constants.F - eSource;
    
end
