function cons = assembleConservationEquation(model, flux, bcflux, source, accum)
% NOTE : scaling with Faraday constant
    
    if nargin < 5
        accum = 0;
    end
        
    op = model.operators;
    F = model.constants.F;
    
    %cons = accum + (op.Div(flux) - bcflux)./model.G.cells.volumes./F - source;
    cons = accum + (op.Div(flux) - bcflux) - source.*model.G.cells.volumes.*F;
end
