function state = setupExchanges(model, state)
    
% We setup the exchange terms between the electrolyte and the electrodes for Li due to chemical reacions:
%   state.Electrolyte.LiSource 
%   state.NegativeElectrode.LiSource 
%   state.PositiveElectrode.LiSource 
%
% We setup the corresponding electron exchange   
%
%   state.Electrolyte.eSource (equal to zero)
%   state.NegativeElectrode.eSource  
%   state.PositiveElectrode.eSource  
%
% We use the following standard abbreviation
%
% elyte : Electrolyte
% ne    : NegativeElectrode
% pe    : PositiveElectrode    
    
    elyte       = model.Electrolyte;
    ne = model.NegativeElectrode;
    pe = model.PositiveElectrode;
    
    elyte_Li_source = zeros(elyte.G.cells.num, 1);
    ne_Li_source    = zeros(ne.G.cells.num, 1);
    pe_Li_source    = zeros(pe.G.cells.num, 1);
    ne_e_source     = zeros(ne.G.cells.num, 1);
    pe_e_source     = zeros(pe.G.cells.num, 1);

    phi = state.Electrolyte.phi;
    if isa(phi, 'ADI')
        adsample = getSampleAD(phi);
        adbackend = model.AutoDiffBackend;
        elyte_Li_source = adbackend.convertToAD(elyte_Li_source, adsample);
        ne_Li_source    = adbackend.convertToAD(ne_Li_source, adsample);
        pe_Li_source    = adbackend.convertToAD(pe_Li_source, adsample);
        ne_e_source     = adbackend.convertToAD(ne_e_source, adsample);
        pe_e_source     = adbackend.convertToAD(pe_e_source, adsample);
    end

    %%%%% Set up chemical source terms %%%%%%%%%%%%%%%%%k

    %%%%% NE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    coupterm = model.getCoupTerm('NegativeElectrode-Electrolyte');
    necells = coupterm.couplingcells(:, 1);
    elytecells = coupterm.couplingcells(:, 2);

    % We compute the reaction rate
    state.NegativeElectrode = ne.updateReactionRate(state.NegativeElectrode);
    ne_R = state.NegativeElectrode.R;

    % elyte NE Li+ source
    elyte_Li_source(elytecells) = ne_R;

    % Active Material NE Li0 source
    ne_Li_source(necells) = - ne_R;

    % Active Material NE current source
    ne_e_source(necells) = + ne_R;

    %%%%% PE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % elyte PE Li+ source
    coupterm = model.getCoupTerm('PositiveElectrode-Electrolyte');
    pecells = coupterm.couplingcells(:, 1);
    elytecells = coupterm.couplingcells(:, 2);

    % calculate rection rate
    state.PositiveElectrode = pe.updateReactionRate(state.PositiveElectrode);
    pe_R = state.PositiveElectrode.R;

    % elyte PE Li+ source
    elyte_Li_source(elytecells) = - pe_R;

    % Active Material PE Li0 source
    pe_Li_source(pecells) = + pe_R;

    % Active Material PE current source
    pe_e_source(pecells) = - pe_R;

    state.Electrolyte.LiSource = elyte_Li_source;
    state.Electrolyte.eSource  = 0;
    state.NegativeElectrode.LiSource = ne_Li_source;
    state.NegativeElectrode.eSource  = ne_e_source;
    state.PositiveElectrode.LiSource = pe_Li_source;
    state.PositiveElectrode.eSource  = pe_e_source;

    state.NegativeCurrentCollector.eSource = 0;
    state.PositiveCurrentCollector.eSource = 0;
    
end