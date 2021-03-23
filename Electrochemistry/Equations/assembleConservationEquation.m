function cons = assembleConservationEquation(model, flux, bcflux, source)
% NOTE : scaling with Faraday constant
    
    op = model.operators;
    F = model.constants.F;
    
    cons = (op.Div(flux) - bcflux)./model.G.cells.volumes./F - source;
    
end
