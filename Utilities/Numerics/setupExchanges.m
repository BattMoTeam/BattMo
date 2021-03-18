function state = setupExchanges(model, state)
    
% We setup the exchange terms between the electrolyte and the electrodes for Li due to chemical reacions:
% Electrolyte_Li_source, NegativeElectrode_Li_source, PositiveElectrode_Li_source.
%
% We setup also the electron production in the electrod: NegativeElectrode_e_source, PositiveElectrode_e_source.
%
                
    Electrolyte = model.Electrolyte;
    NegativeElectrode    = model.NegativeElectrode;
    PositiveElectrode    = model.PositiveElectrode;
    
    Electrolyte_Li_source = zeros(Electrolyte.G.cells.num, 1);
    NegativeElectrode_Li_source    = zeros(NegativeElectrode.G.cells.num, 1);
    PositiveElectrode_Li_source    = zeros(PositiveElectrode.G.cells.num, 1);
    NegativeElectrode_e_source     = zeros(NegativeElectrode.G.cells.num, 1);
    PositiveElectrode_e_source     = zeros(PositiveElectrode.G.cells.num, 1);

    phi = state.Electrolyte.phi;
    if isa(phi, 'ADI')
        adsample = getSampleAD(phi);
        adbackend = model.AutoDiffBackend;
        Electrolyte_Li_source = adbackend.convertToAD(Electrolyte_Li_source, adsample);
        NegativeElectrode_Li_source    = adbackend.convertToAD(NegativeElectrode_Li_source, adsample);
        PositiveElectrode_Li_source    = adbackend.convertToAD(PositiveElectrode_Li_source, adsample);
        NegativeElectrode_e_source     = adbackend.convertToAD(NegativeElectrode_e_source, adsample);
        PositiveElectrode_e_source     = adbackend.convertToAD(PositiveElectrode_e_source, adsample);
    end

    %%%%% Set up chemical source terms %%%%%%%%%%%%%%%%%k

    %%%%% NE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    coupterm = model.getCoupTerm('NegativeElectrode-Electrolyte');
    NegativeElectrodecells = coupterm.couplingcells(:, 1);
    Electrolytecells = coupterm.couplingcells(:, 2);

    % We compute the reaction rate
    state.NegativeElectrode = NegativeElectrode.updateReactionRate(state.NegativeElectrode);
    NegativeElectrode_R = state.NegativeElectrode.R;

    % Electrolyte NE Li+ source
    Electrolyte_Li_source(Electrolytecells) = NegativeElectrode_R;

    % Active Material NE Li0 source
    NegativeElectrode_Li_source(NegativeElectrodecells) = - NegativeElectrode_R;

    % Active Material NE current source
    NegativeElectrode_e_source(NegativeElectrodecells) = + NegativeElectrode_R;

    %%%%% PE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Electrolyte PE Li+ source
    coupterm = model.getCoupTerm('PositiveElectrode-Electrolyte');
    PositiveElectrodecells = coupterm.couplingcells(:, 1);
    Electrolytecells = coupterm.couplingcells(:, 2);

    % calculate rection rate
    state.PositiveElectrode = PositiveElectrode.updateReactionRate(state.PositiveElectrode);
    PositiveElectrode_R = state.PositiveElectrode.R;

    % Electrolyte PE Li+ source
    Electrolyte_Li_source(Electrolytecells) = - PositiveElectrode_R;

    % Active Material PE Li0 source
    PositiveElectrode_Li_source(PositiveElectrodecells) = + PositiveElectrode_R;

    % Active Material PE current source
    PositiveElectrode_e_source(PositiveElectrodecells) = - PositiveElectrode_R;

    state.Electrolyte.LiSource =  Electrolyte_Li_source;
    state.NegativeElectrode.LiSource    =  NegativeElectrode_Li_source;
    state.NegativeElectrode.eSource     =  NegativeElectrode_e_source;
    state.PositiveElectrode.LiSource    =  PositiveElectrode_Li_source;
    state.PositiveElectrode.eSource     =  PositiveElectrode_e_source;
    
end