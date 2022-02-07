function [mass, masses] = computeCellMass(model)

    
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
        

        masses.(elde).(eac).val = sum(rho.*vols.*poro);
        
        mass = mass + masses.(elde).(eac).val;
        
        rho  = model.(elde).(cc).density;
        vols = model.(elde).(cc).G.cells.volumes;
        
        masses.(elde).(cc).val = sum(rho.*vols);
        mass = mass + masses.(elde).(cc).val;
        
    end
    
    rho  = model.(elyte).density;
    vols = model.(elyte).G.cells.volumes;
    poro = model.(elyte).volumeFraction;
    
    masses.(elyte).val = sum(rho.*vols.*poro);
    mass = mass + masses.(elyte).val;
    
    rho  = model.(elyte).(sep).density;
    vols = model.(elyte).(sep).G.cells.volumes;
    poro = model.(elyte).(sep).volumeFraction;
    
    masses.(elyte).(sep).val = sum(rho.*vols.*poro);
    mass = mass + masses.(elyte).(sep).val;
    
end
