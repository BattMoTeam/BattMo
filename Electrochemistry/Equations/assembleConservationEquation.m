function cons = assembleConservationEquation(model, flux, bcflux, source)
    
    op = model.operators;
    cons = (op.Div(flux) - bcflux)./ model.G.cells.volumes./model.constants.F - source;
    
end
