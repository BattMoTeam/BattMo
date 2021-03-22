function flux = assembleFlux(model, potential, fluxCoefficient)
    
    op = model.operators;
    flux = - op.harmFace(fluxCoefficient).*op.Grad(phi); 
    
end
