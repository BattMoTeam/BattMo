function  jHeatOhmSource = updateOhmSourceFunc(model, j, effectiveElectricalConductivity)
%
% reference : 
% @misc{Latz_2016,
% 	doi = {10.1515/nano.bjneah.6.102},
% 	url = {https://doi.org/10.1515%2Fnano.bjneah.6.102},
% 	year = 2016,
% 	month = {jul},
% 	publisher = {De Gruyter},
% 	author = {Arnulf Latz and Jochen Zausch},
% 	title = {Multiscale modeling of lithium ion batteries: thermal aspects}
% }
    
    op = model.operators.cellFluxOp;
    
    if isprop(model, 'volumeFraction')
        volfrac = model.volumeFraction;
    else
        volfrac = 1;
    end
    vols = model.G.cells.volumes;

    j = op.P*j;
    jsq = j.^2;
    jsq = op.S*jsq;
    
    
    jHeatOhmSource = (volfrac.*vols).*jsq./effectiveElectricalConductivity;
    
    
end