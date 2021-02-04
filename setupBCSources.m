function state = setupBCSources(model, state)
    
    ne    = model.getAssocModel('ne');
    pe    = model.getAssocModel('pe');
    ccne  = model.getAssocModel('ccne');
    ccpe  = model.getAssocModel('ccpe');
    
    ne_am = ne.getAssocModel('am');
    pe_am = pe.getAssocModel('am');
    
    ne_phi   = state.ne.am.phi;
    pe_phi   = state.pe.am.phi;
    ccne_phi = state.ccne.phi;
    ccpe_phi = state.ccpe.phi;
    ccpe_E   = state.ccpe.E;
    
    
    %% We assume state.ne.am.D and state.ne.am.D  have been updated

    D = state.ne.am.D;
    ne_Deff = D.*ne_am.eps.^1.5;

    D = state.pe.am.D;
    pe_Deff = D.*pe_am.eps.^1.5;

    ne_sigmaeff = ne.sigmaeff;
    pe_sigmaeff = pe.sigmaeff;
    ccne_sigmaeff = ccne.sigmaeff;
    ccpe_sigmaeff = ccpe.sigmaeff;

    %% We setup the current transfers between ccne collector and ne material with ne_j_bcsource and ccne_j_bcsource

    ne_j_bcsource = ne_phi*0.0; %NB hack to initialize zero ad
    ccne_j_bcsource = ccne_phi*0.0; %NB hack to initialize zero ad

    coupterm = model.getCoupTerm('ccne-ne');
    face_ccne = coupterm.couplingfaces(:, 1);
    face_ne = coupterm.couplingfaces(:, 2);
    [tne, bccell_ne] = ne.operators.harmFaceBC(ne_sigmaeff, face_ne);
    [tccne, bccell_ccne] = ccne.operators.harmFaceBC(ccne_sigmaeff, face_ccne);

    bcphi_ne = ne_phi(bccell_ne);
    bcphi_ccne = ccne_phi(bccell_ccne);

    trans = 1./(1./tne + 1./tccne);
    crosscurrent = trans.*(bcphi_ccne - bcphi_ne);
    ne_j_bcsource(bccell_ne) = crosscurrent;
    ccne_j_bcsource(bccell_ccne) = -crosscurrent;

    %% We impose the boundary condition at chosen boundary cells of the ne current collector by updating ccne_j_bcsource

    coupterm = model.getCoupTerm('bc-ccne');
    faces = coupterm.couplingfaces;
    bcval = zeros(numel(faces), 1);
    [tccne, cells] = ccne.operators.harmFaceBC(ccne_sigmaeff, faces);
    ccne_j_bcsource(cells) = ccne_j_bcsource(cells) + tccne.*(bcval - ccne_phi(cells));

    %% We setup the current transfers between ccpe collector and pe material with pe_j_bcsource and ccpe_j_bcsource

    pe_j_bcsource   = pe_phi*0.0; %NB hack to initialize zero ad
    ccpe_j_bcsource = ccpe_phi*0.0; %NB hack to initialize zero ad

    coupterm = model.getCoupTerm('ccpe-pe');
    face_ccpe = coupterm.couplingfaces(:, 1);
    face_pe = coupterm.couplingfaces(:, 2);
    [tpe, bccell_pe] = pe.operators.harmFaceBC(pe_sigmaeff, face_pe);
    [tccpe, bccell_ccpe] = ccpe.operators.harmFaceBC(ccpe_sigmaeff, face_ccpe);
    bcphi_pe = pe_phi(bccell_pe);
    bcphi_ccpe = ccpe_phi(bccell_ccpe);

    trans = 1./(1./tpe + 1./tccpe);
    crosscurrent = trans.*(bcphi_ccpe - bcphi_pe);
    pe_j_bcsource(bccell_pe) = crosscurrent;
    ccpe_j_bcsource(bccell_ccpe) = -crosscurrent;

    %% We impose the boundary condition at chosen boundary cells of the anode current collector by updating ccpe_j_bcsource

    coupterm = model.getCoupTerm('bc-ccpe');
    faces = coupterm.couplingfaces;
    bcval = ccpe_E;
    [tccpe, cells] = ccpe.operators.harmFaceBC(ccpe_sigmaeff, faces);
    ccpe_j_bcsource(cells) = ccpe_j_bcsource(cells) + tccpe.*(bcval - ccpe_phi(cells));
    
    state.ne.jBcSource   =  ne_j_bcsource;
    state.pe.jBcSource   =  pe_j_bcsource;
    state.ccne.jBcSource =  ccne_j_bcsource;
    state.ccpe.jBcSource =  ccpe_j_bcsource;
    
end