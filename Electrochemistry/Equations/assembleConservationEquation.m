function cons = assembleConservationEquation(model, flux, bcsource, source, accum)
    
    if nargin < 5
        accum = 0;
    end
        
    op = model.operators;
    
    % cons = accum + (op.Div(flux) - bcsource) - source;
    accum = accum - bcsource - source; 
    cons = op.AccDiv(accum, flux);
    %cons = cons - bcsource - source;
    
end
