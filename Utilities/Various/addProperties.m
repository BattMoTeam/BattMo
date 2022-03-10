function state = addProperties(model, state)

    nc = model.G.cells.num;
    %state.SOC = model.SOC*ones(nc, 1);

    % Shorthands used in this function
    battery = model;
    ne      = 'NegativeElectrode';
    pe      = 'PositiveElectrode';
    am     = 'ActiveMaterial';
    cc      = 'CurrentCollector';
    elyte   = 'Electrolyte';
    am      = 'ActiveMaterial';
    thermal = 'ThermalModel';

    electrodes = {ne, pe};
    electrodecomponents = {am, cc};

    %% Synchronization across components

    % temperature
    state = battery.updateTemperature(state);

    state.(elyte) = battery.(elyte).updateConcentrations(state.(elyte));

    for ind = 1 : numel(electrodes)
        elde = electrodes{ind};
        % potential and concentration between active material and electode active component
        state.(elde).(am) = battery.(elde).(am).updatePhi(state.(elde).(am));
        if(model.use_solid_diffusion)
            state.(elde).(am) = battery.(elde).(am).updateChargeCarrier(state.(elde).(am));
        else
            state.(elde).(am).c = state.(elde).(am).(am).cElectrode;
            state.(elde).(am) = battery.(elde).(am).updateChargeCarrier(state.(elde).(am));                    
        end              
    end


    %% Update Electrolyte -> Electrodes coupling 

    state = battery.updateElectrodeCoupling(state); 

    %% Update reaction rates in both electrodes

    for ind = 1 : numel(electrodes)
        elde = electrodes{ind};
        state.(elde).(am).(am) = battery.(elde).(am).(am).updateReactionRateCoefficient(state.(elde).(am).(am));
        state.(elde).(am).(am) = battery.(elde).(am).(am).updateDiffusionCoefficient(state.(elde).(am).(am));
        state.(elde).(am).(am) = battery.(elde).(am).(am).updateOCP(state.(elde).(am).(am));
        state.(elde).(am).(am) = battery.(elde).(am).(am).updateReactionRate(state.(elde).(am).(am));
    end

    %% Update Electrodes -> Electrolyte  coupling

    state = battery.updateElectrolyteCoupling(state);

    %% Update coupling within electrodes and external coupling

    state.(ne) = battery.(ne).updateCoupling(state.(ne));
    state.(pe) = battery.(pe).updateCoupling(state.(pe));

    state.(ne).(am) = battery.(ne).(am).updatejBcSource(state.(ne).(am));
    state.(pe).(am) = battery.(pe).(am).updatejBcSource(state.(pe).(am));

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
        state.(elde).(am) = battery.(elde).(am).updateIonAndCurrentSource(state.(elde).(am));
        state.(elde).(am) = battery.(elde).(am).updateCurrent(state.(elde).(am));
        state.(elde).(am) = battery.(elde).(am).updateChargeConservation(state.(elde).(am));
    end

    %% elyte mass conservation

    state.(elyte) = battery.(elyte).updateDiffusionCoefficient(state.(elyte));
    state.(elyte) = battery.(elyte).updateMassFlux(state.(elyte));


    for ind = 1 : numel(electrodes)
        elde = electrodes{ind};
        
        %% Electrodes mass conservation
        state.(elde).(am) = battery.(elde).(am).updateMassFlux(state.(elde).(am));
        
        %% Electrodes charge conservation - current collector part
        state.(elde).(cc) = battery.(elde).(cc).updateCurrent(state.(elde).(cc));
        state.(elde).(cc) = battery.(elde).(cc).updateChargeConservation(state.(elde).(cc));

    end

    %% update solid diffustion equations
    for ind = 1 : numel(electrodes)
        elde = electrodes{ind};
        state.(elde).(am).(am) = battery.(elde).(am).(am).assembleSolidDiffusionEquation(state.(elde).(am).(am));
    end

    %% update Face fluxes
    for ind = 1 : numel(electrodes)
        elde = electrodes{ind};
        state.(elde).(am) = battery.(elde).(am).updateFaceCurrent(state.(elde).(am));
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


%{
Copyright 2009-2021 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
