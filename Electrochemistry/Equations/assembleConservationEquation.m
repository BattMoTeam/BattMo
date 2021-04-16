function cons = assembleConservationEquation(model, flux, bcflux, source, accum)
    
    if nargin < 5
        accum = 0;
    end
        
    op = model.operators;
    F = model.constants.F;
    
    cons = accum + (op.Div(flux) - bcflux) - source;
    
end
