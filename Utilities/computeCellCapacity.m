function c = computeCellCapacity(model)

%
%
% SYNOPSIS:
%   function c = computeCellCapacity(model, state)
%
% DESCRIPTION: computes the cell usable capacity
%
% PARAMETERS:
%   model - battery model
%
% RETURNS:
%   c - capacity
%
% EXAMPLE:
%
% SEE ALSO:
%

% We handle the battery and bare battery structure (without current collector)
    
    ne  = 'NegativeElectrode';
    pe  = 'PositiveElectrode';
    eac = 'ElectrodeActiveComponent';
    am  = 'ActiveMaterial';
    
    eldes = {ne, pe};
    
    isBare = false;
    if strcmp(class(model), 'BareBattery')
        isBare = true;
    end
    
    for ind = 1 : numel(eldes)
        
        elde = eldes{ind};
        
        if isBare
            ammodel = model.(elde).(am);
        else
            ammodel = model.(elde).(eac).(am);
        end
        
        n = ammodel.n;
        F = ammodel.constants.F;
        G = ammodel.G;
        c_max    = ammodel.Li.cmax;
        theta0   = ammodel.theta0;
        theta100 = ammodel.theta100;
        volume_fraction = ammodel.volumeFraction;
        
        volume_electrode = sum(G.cells.volumes);
        
        c_usable = abs(theta100 - theta0) * c_max;
        mol_usable = c_usable * volume_fraction * volume_electrode;
        
        cap_usable(ind) = mol_usable * n * F;
        
    end
    
    c = min(cap_usable); 
    
end
