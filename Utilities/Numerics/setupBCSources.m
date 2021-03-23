function state = setupBCSources(model, state)
% Abbreviations used in this function:
% elyte : Electrolyte
% ne    : NegativeElectrode
% pe    : PositiveElectrode
% ncc   : NegativeCurrentCollector
% pcc   : PositiveCurrentCollector
    
    elyte = model.Electrolyte;
    ne  = model.NegativeElectrode;
    pe  = model.PositiveElectrode;
    ncc = model.NegativeCurrentCollector;
    pcc = model.PositiveCurrentCollector;
    
    ne_am = ne.ActiveMaterial;
    pe_am = pe.ActiveMaterial;
    
    ne_phi  = state.NegativeElectrode.ActiveMaterial.phi;
    pe_phi  = state.PositiveElectrode.ActiveMaterial.phi;
    ncc_phi = state.NegativeCurrentCollector.phi;
    pcc_phi = state.PositiveCurrentCollector.phi;
    pcc_E   = state.PositiveCurrentCollector.E;
    
    %% We assume state.NegativeElectrode.ActiveMaterial.D and state.NegativeElectrode.ActiveMaterial.D  have been updated

    D = state.NegativeElectrode.ActiveMaterial.D;
    ne_Deff = D.*ne_am.volumeFraction.^1.5;

    D = state.PositiveElectrode.ActiveMaterial.D;
    pe_Deff = D.*pe_am.volumeFraction.^1.5;

    ne_sigmaeff  = ne.EffectiveElectronicConductivity;
    pe_sigmaeff  = pe.EffectiveElectronicConductivity;
    ncc_sigmaeff = ncc.EffectiveElectronicConductivity;
    pcc_sigmaeff = pcc.EffectiveElectronicConductivity;

    %% We setup the current transfers between NegativeCurrentCollector collector and ne material with ne_j_bcsource and ncc_j_bcsource

    ne_j_bcsource  = ne_phi*0.0; %NB hack to initialize zero ad
    ncc_j_bcsource = ncc_phi*0.0; %NB hack to initialize zero ad

    coupterm = model.getCoupTerm('NegativeCurrentCollector-NegativeElectrode');
    face_ncc = coupterm.couplingfaces(:, 1);
    face_ne = coupterm.couplingfaces(:, 2);
    [tne, bccell_ne] = ne.operators.harmFaceBC(ne_sigmaeff, face_ne);
    [tncc, bccell_ncc] = ncc.operators.harmFaceBC(ncc_sigmaeff, face_ncc);

    bcphi_ne = ne_phi(bccell_ne);
    bcphi_ncc = ncc_phi(bccell_ncc);

    trans = 1./(1./tne + 1./tncc);
    crosscurrent = trans.*(bcphi_ncc - bcphi_ne);
    ne_j_bcsource(bccell_ne) = crosscurrent;
    ncc_j_bcsource(bccell_ncc) = -crosscurrent;

    %% We impose the boundary condition at chosen boundary cells of the NegativeElectrode current collector by updating ncc_j_bcsource

    coupterm = model.getCoupTerm('bc-NegativeCurrentCollector');
    faces = coupterm.couplingfaces;
    bcval = zeros(numel(faces), 1);
    [tncc, cells] = ncc.operators.harmFaceBC(ncc_sigmaeff, faces);
    ncc_j_bcsource(cells) = ncc_j_bcsource(cells) + tncc.*(bcval - ncc_phi(cells));

    %% We setup the current transfers between PositiveCurrentCollector collector and PositiveElectrode material with pe_j_bcsource and pcc_j_bcsource

    pe_j_bcsource   = pe_phi*0.0; %NB hack to initialize zero ad
    pcc_j_bcsource = pcc_phi*0.0; %NB hack to initialize zero ad

    coupterm = model.getCoupTerm('PositiveCurrentCollector-PositiveElectrode');
    face_pcc = coupterm.couplingfaces(:, 1);
    face_pe = coupterm.couplingfaces(:, 2);
    [tpe, bccell_pe] = pe.operators.harmFaceBC(pe_sigmaeff, face_pe);
    [tpcc, bccell_pcc] = pcc.operators.harmFaceBC(pcc_sigmaeff, face_pcc);
    bcphi_pe = pe_phi(bccell_pe);
    bcphi_pcc = pcc_phi(bccell_pcc);

    trans = 1./(1./tpe + 1./tpcc);
    crosscurrent = trans.*(bcphi_pcc - bcphi_pe);
    pe_j_bcsource(bccell_pe) = crosscurrent;
    pcc_j_bcsource(bccell_pcc) = -crosscurrent;

    %% We impose the boundary condition at chosen boundary cells of the anode current collector by updating PositiveCurrentCollector_j_bcsource

    coupterm = model.getCoupTerm('bc-PositiveCurrentCollector');
    faces = coupterm.couplingfaces;
    bcval = pcc_E;
    [tpcc, cells] = pcc.operators.harmFaceBC(pcc_sigmaeff, faces);
    pcc_j_bcsource(cells) = pcc_j_bcsource(cells) + tpcc.*(bcval - pcc_phi(cells));
    
    state.Electrolyte.jBcSource = 0;
    state.NegativeElectrode.jBcSource =  ne_j_bcsource;
    state.PositiveElectrode.jBcSource =  pe_j_bcsource;
    state.NegativeCurrentCollector.jBcSource =  ncc_j_bcsource;
    state.PositiveCurrentCollector.jBcSource =  pcc_j_bcsource;
   
end