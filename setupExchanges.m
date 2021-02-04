function state = setupExchanges(model, state)
    
% We setup the exchange terms between the electrolyte and the electrodes for Li due to chemical reacions:
% elyte_Li_source, ne_Li_source, pe_Li_source.
%
% We setup also the electron production in the electrod: ne_e_source, pe_e_source.
%
                
    elyte = model.elyte;
    ne    = model.ne;
    pe    = model.pe;
    
    elyte_Li_source = zeros(elyte.G.cells.num, 1);
    ne_Li_source    = zeros(ne.G.cells.num, 1);
    pe_Li_source    = zeros(pe.G.cells.num, 1);
    ne_e_source     = zeros(ne.G.cells.num, 1);
    pe_e_source     = zeros(pe.G.cells.num, 1);

    phi = state.elyte.phi;
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

    coupterm = model.getCoupTerm('ne-elyte');
    necells = coupterm.couplingcells(:, 1);
    elytecells = coupterm.couplingcells(:, 2);

    % We compute the reaction rate
    state.ne = ne.updateReactionRate(state.ne);
    ne_R = state.ne.R;

    % Electrolyte NE Li+ source
    elyte_Li_source(elytecells) = ne_R;

    % Active Material NE Li0 source
    ne_Li_source(necells) = - ne_R;

    % Active Material NE current source
    ne_e_source(necells) = + ne_R;

    %%%%% PE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Electrolyte PE Li+ source
    coupterm = model.getCoupTerm('pe-elyte');
    pecells = coupterm.couplingcells(:, 1);
    elytecells = coupterm.couplingcells(:, 2);

    % calculate rection rate
    state.pe = pe.updateReactionRate(state.pe);
    pe_R = state.pe.R;

    % Electrolyte PE Li+ source
    elyte_Li_source(elytecells) = - pe_R;

    % Active Material PE Li0 source
    pe_Li_source(pecells) = + pe_R;

    % Active Material PE current source
    pe_e_source(pecells) = - pe_R;

    state.elyte.LiSource =  elyte_Li_source;
    state.ne.LiSource    =  ne_Li_source;
    state.ne.eSource     =  ne_e_source;
    state.pe.LiSource    =  pe_Li_source;
    state.pe.eSource     =  pe_e_source;
    
end