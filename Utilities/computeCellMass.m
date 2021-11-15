function mass = computeCellMass(model)

    
    elyte = 'Electrolyte';
    sep   = 'Separator';
    ne    = 'NegativeElectrode';
    pe    = 'PositiveElectrode';
    eac   = 'ElectrodeActiveComponent';
    am    = 'ActiveMaterial';
    cc    = 'CurrentCollector';
    
    eldes = {ne, pe};
    
    mass = 0; % total mass
    
    % We first compute the mass of the electrodes
    for ind = 1 : numel(eldes)
        
        elde = eldes{ind};
        
        rho  = model.(elde).(eac).(am).density;
        vols = model.(elde).(eac).G.cells.volumes;
        poro = model.(elde).(eac).volumeFraction;
        
        mass = mass + sum(rho.*vols.*poro);
        
        rho  = model.(elde).(cc).density;
        vols = model.(elde).(cc).G.cells.volumes;
        
        mass = mass + sum(rho.*vols);
        
    end
    
    rho  = model.(elyte).density;
    vols = model.(elyte).G.cells.volumes;
    poro = model.(elyte).volumeFraction;
    
    mass = mass + sum(rho.*vols.*poro);
    
    rho  = model.(elyte).(sep).density;
    vols = model.(elyte).(sep).G.cells.volumes;
    poro = model.(elyte).(sep).volumeFraction;
    
    mass = mass + sum(rho.*vols.*poro);
    
    
end
