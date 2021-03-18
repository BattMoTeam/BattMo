function state = setupBCSources(model, state)
    
    NegativeElectrode    = model.NegativeElectrode;
    PositiveElectrode    = model.PositiveElectrode;
    NegativeCurrentCollector  = model.NegativeCurrentCollector;
    PositiveCurrentCollector  = model.PositiveCurrentCollector;
    
    NegativeElectrode_ActiveMaterial = NegativeElectrode.ActiveMaterial;
    PositiveElectrode_ActiveMaterial = PositiveElectrode.ActiveMaterial;
    
    NegativeElectrode_phi   = state.NegativeElectrode.ActiveMaterial.phi;
    PositiveElectrode_phi   = state.PositiveElectrode.ActiveMaterial.phi;
    NegativeCurrentCollector_phi = state.NegativeCurrentCollector.phi;
    PositiveCurrentCollector_phi = state.PositiveCurrentCollector.phi;
    PositiveCurrentCollector_E   = state.PositiveCurrentCollector.E;
    
    
    %% We assume state.NegativeElectrode.ActiveMaterial.D and state.NegativeElectrode.ActiveMaterial.D  have been updated

    D = state.NegativeElectrode.ActiveMaterial.D;
    NegativeElectrode_Deff = D.*NegativeElectrode_ActiveMaterial.volumeFraction.^1.5;

    D = state.PositiveElectrode.ActiveMaterial.D;
    PositiveElectrode_Deff = D.*PositiveElectrode_ActiveMaterial.volumeFraction.^1.5;

    NegativeElectrode_sigmaeff = NegativeElectrode.effectiveElectronicConductivity;
    PositiveElectrode_sigmaeff = PositiveElectrode.effectiveElectronicConductivity;
    NegativeCurrentCollector_sigmaeff = NegativeCurrentCollector.effectiveElectronicConductivity;
    PositiveCurrentCollector_sigmaeff = PositiveCurrentCollector.effectiveElectronicConductivity;

    %% We setup the current transfers between NegativeCurrentCollector collector and NegativeElectrode material with NegativeElectrode_j_bcsource and NegativeCurrentCollector_j_bcsource

    NegativeElectrode_j_bcsource = NegativeElectrode_phi*0.0; %NB hack to initialize zero ad
    NegativeCurrentCollector_j_bcsource = NegativeCurrentCollector_phi*0.0; %NB hack to initialize zero ad

    coupterm = model.getCoupTerm('NegativeCurrentCollector-NegativeElectrode');
    face_NegativeCurrentCollector = coupterm.couplingfaces(:, 1);
    face_NegativeElectrode = coupterm.couplingfaces(:, 2);
    [tNegativeElectrode, bccell_NegativeElectrode] = NegativeElectrode.operators.harmFaceBC(NegativeElectrode_sigmaeff, face_NegativeElectrode);
    [tNegativeCurrentCollector, bccell_NegativeCurrentCollector] = NegativeCurrentCollector.operators.harmFaceBC(NegativeCurrentCollector_sigmaeff, face_NegativeCurrentCollector);

    bcphi_NegativeElectrode = NegativeElectrode_phi(bccell_NegativeElectrode);
    bcphi_NegativeCurrentCollector = NegativeCurrentCollector_phi(bccell_NegativeCurrentCollector);

    trans = 1./(1./tNegativeElectrode + 1./tNegativeCurrentCollector);
    crosscurrent = trans.*(bcphi_NegativeCurrentCollector - bcphi_NegativeElectrode);
    NegativeElectrode_j_bcsource(bccell_NegativeElectrode) = crosscurrent;
    NegativeCurrentCollector_j_bcsource(bccell_NegativeCurrentCollector) = -crosscurrent;

    %% We impose the boundary condition at chosen boundary cells of the NegativeElectrode current collector by updating NegativeCurrentCollector_j_bcsource

    coupterm = model.getCoupTerm('bc-NegativeCurrentCollector');
    faces = coupterm.couplingfaces;
    bcval = zeros(numel(faces), 1);
    [tNegativeCurrentCollector, cells] = NegativeCurrentCollector.operators.harmFaceBC(NegativeCurrentCollector_sigmaeff, faces);
    NegativeCurrentCollector_j_bcsource(cells) = NegativeCurrentCollector_j_bcsource(cells) + tNegativeCurrentCollector.*(bcval - NegativeCurrentCollector_phi(cells));

    %% We setup the current transfers between PositiveCurrentCollector collector and PositiveElectrode material with PositiveElectrode_j_bcsource and PositiveCurrentCollector_j_bcsource

    PositiveElectrode_j_bcsource   = PositiveElectrode_phi*0.0; %NB hack to initialize zero ad
    PositiveCurrentCollector_j_bcsource = PositiveCurrentCollector_phi*0.0; %NB hack to initialize zero ad

    coupterm = model.getCoupTerm('PositiveCurrentCollector-PositiveElectrode');
    face_PositiveCurrentCollector = coupterm.couplingfaces(:, 1);
    face_PositiveElectrode = coupterm.couplingfaces(:, 2);
    [tPositiveElectrode, bccell_PositiveElectrode] = PositiveElectrode.operators.harmFaceBC(PositiveElectrode_sigmaeff, face_PositiveElectrode);
    [tPositiveCurrentCollector, bccell_PositiveCurrentCollector] = PositiveCurrentCollector.operators.harmFaceBC(PositiveCurrentCollector_sigmaeff, face_PositiveCurrentCollector);
    bcphi_PositiveElectrode = PositiveElectrode_phi(bccell_PositiveElectrode);
    bcphi_PositiveCurrentCollector = PositiveCurrentCollector_phi(bccell_PositiveCurrentCollector);

    trans = 1./(1./tPositiveElectrode + 1./tPositiveCurrentCollector);
    crosscurrent = trans.*(bcphi_PositiveCurrentCollector - bcphi_PositiveElectrode);
    PositiveElectrode_j_bcsource(bccell_PositiveElectrode) = crosscurrent;
    PositiveCurrentCollector_j_bcsource(bccell_PositiveCurrentCollector) = -crosscurrent;

    %% We impose the boundary condition at chosen boundary cells of the anode current collector by updating PositiveCurrentCollector_j_bcsource

    coupterm = model.getCoupTerm('bc-PositiveCurrentCollector');
    faces = coupterm.couplingfaces;
    bcval = PositiveCurrentCollector_E;
    [tPositiveCurrentCollector, cells] = PositiveCurrentCollector.operators.harmFaceBC(PositiveCurrentCollector_sigmaeff, faces);
    PositiveCurrentCollector_j_bcsource(cells) = PositiveCurrentCollector_j_bcsource(cells) + tPositiveCurrentCollector.*(bcval - PositiveCurrentCollector_phi(cells));
    
    state.NegativeElectrode.jBcSource   =  NegativeElectrode_j_bcsource;
    state.PositiveElectrode.jBcSource   =  PositiveElectrode_j_bcsource;
    state.NegativeCurrentCollector.jBcSource =  NegativeCurrentCollector_j_bcsource;
    state.PositiveCurrentCollector.jBcSource =  PositiveCurrentCollector_j_bcsource;
    
end