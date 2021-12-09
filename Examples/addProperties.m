function state = addProperties(model, state)

    nc = model.G.cells.num;
    %state.SOC = model.SOC*ones(nc, 1);

    % Shortcuts used in this function
    battery = model;
    ne      = 'NegativeElectrode';
    pe      = 'PositiveElectrode';
    eac     = 'ElectrodeActiveComponent';
    cc      = 'CurrentCollector';
    elyte   = 'Electrolyte';
    am      = 'ActiveMaterial';
    thermal = 'ThermalModel';

    electrodes = {ne, pe};
    electrodecomponents = {eac, cc};

    %% Synchronization across components

    % temperature
    state = battery.updateTemperature(state);

    state.(elyte) = battery.(elyte).updateConcentrations(state.(elyte));

    for ind = 1 : numel(electrodes)
        elde = electrodes{ind};
        % potential and concentration between active material and electode active component
        state.(elde).(eac) = battery.(elde).(eac).updatePhi(state.(elde).(eac));
        if(model.use_solid_diffusion)
            state.(elde).(eac) = battery.(elde).(eac).updateChargeCarrier(state.(elde).(eac));
        else
            state.(elde).(eac).c = state.(elde).(eac).(am).cElectrode;
            state.(elde).(eac) = battery.(elde).(eac).updateChargeCarrier(state.(elde).(eac));                    
        end              
    end


    %% Update Electrolyte -> Electrodes coupling 

    state = battery.updateElectrodeCoupling(state); 

    %% Update reaction rates in both electrodes

    for ind = 1 : numel(electrodes)
        elde = electrodes{ind};
        state.(elde).(eac).(am) = battery.(elde).(eac).(am).updateReactionRateCoefficient(state.(elde).(eac).(am));
        state.(elde).(eac).(am) = battery.(elde).(eac).(am).updateDiffusionCoefficient(state.(elde).(eac).(am));
        state.(elde).(eac).(am) = battery.(elde).(eac).(am).updateOCP(state.(elde).(eac).(am));
        state.(elde).(eac).(am) = battery.(elde).(eac).(am).updateReactionRate(state.(elde).(eac).(am));
    end

    %% Update Electrodes -> Electrolyte  coupling

    state = battery.updateElectrolyteCoupling(state);

    %% Update coupling within electrodes and external coupling

    state.(ne) = battery.(ne).updateCoupling(state.(ne));
    state.(pe) = battery.(pe).updateCoupling(state.(pe));

    state.(ne).(eac) = battery.(ne).(eac).updatejBcSource(state.(ne).(eac));
    state.(pe).(eac) = battery.(pe).(eac).updatejBcSource(state.(pe).(eac));

    state = model.setupExternalCouplingNegativeElectrode(state);
    state = model.setupExternalCouplingPositiveElectrode(state);

    state.(ne).(cc) = battery.(ne).(cc).updatejBcSource(state.(ne).(cc));
    state.(pe).(cc) = battery.(pe).(cc).updatejBcSource(state.(pe).(cc));

    %% elyte charge conservation

    state.(elyte) = battery.(elyte).updateCurrentBcSource(state.(elyte));
    state.(elyte) = battery.(elyte).updateConductivity(state.(elyte));
    state.(elyte) = battery.(elyte).updateChemicalCurrent(state.(elyte));
    state.(elyte) = battery.(elyte).updateCurrent(state.(elyte));
    state.(elyte) = battery.(elyte).updateChargeConservation(state.(elyte));

    %% Electrodes charge conservation - Active material part

    for ind = 1 : numel(electrodes)
        elde = electrodes{ind};
        state.(elde).(eac) = battery.(elde).(eac).updateIonAndCurrentSource(state.(elde).(eac));
        state.(elde).(eac) = battery.(elde).(eac).updateCurrent(state.(elde).(eac));
        state.(elde).(eac) = battery.(elde).(eac).updateChargeConservation(state.(elde).(eac));
    end

    %% elyte mass conservation

    state.(elyte) = battery.(elyte).updateDiffusionCoefficient(state.(elyte));
    state.(elyte) = battery.(elyte).updateMassFlux(state.(elyte));


    for ind = 1 : numel(electrodes)
        elde = electrodes{ind};
        
        %% Electrodes mass conservation
        state.(elde).(eac) = battery.(elde).(eac).updateMassFlux(state.(elde).(eac));
        
        %% Electrodes charge conservation - current collector part
        state.(elde).(cc) = battery.(elde).(cc).updateCurrent(state.(elde).(cc));
        state.(elde).(cc) = battery.(elde).(cc).updateChargeConservation(state.(elde).(cc));

    end

    %% update solid diffustion equations
    for ind = 1 : numel(electrodes)
        elde = electrodes{ind};
        state.(elde).(eac).(am) = battery.(elde).(eac).(am).assembleSolidDiffusionEquation(state.(elde).(eac).(am));
    end

    %% update Face fluxes
    for ind = 1 : numel(electrodes)
        elde = electrodes{ind};
        state.(elde).(eac) = battery.(elde).(eac).updateFaceCurrent(state.(elde).(eac));
        state.(elde).(cc) = battery.(elde).(cc).updateFaceCurrent(state.(elde).(cc));
    end
    state.(elyte) = battery.(elyte).updateFaceCurrent(state.(elyte));

    %% update Thermal source term from electrical resistance

    state = battery.updateThermalOhmicSourceTerms(state);
    state = battery.updateThermalChemicalSourceTerms(state);
    state = battery.updateThermalReactionSourceTerms(state);

    state.(thermal) = battery.(thermal).updateHeatSourceTerm(state.(thermal));
    state.(thermal) = battery.(thermal).updateThermalBoundarySourceTerms(state.(thermal));
    
end