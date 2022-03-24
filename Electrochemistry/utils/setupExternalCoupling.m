function [jExternal, jFaceExternal] = setupExternalCoupling(model, phi, phiExternal)
    
    coupterm = model.externalCouplingTerm;
    
    jExternal = phi*0.0; %NB hack to initialize zero ad
    
    sigmaeff = model.EffectiveElectricalConductivity;
    faces = coupterm.couplingfaces;
    bcval = phiExternal;
    [t, cells] = model.operators.harmFaceBC(sigmaeff, faces);
    current = t.*(bcval - phi(cells));
    %jExternal(cells) = jExternal(cells) + current;
    jExternal = subsetPlus(jExternal,current,cells);
    G = model.G;
    nf = G.faces.num;
    sgn = model.operators.sgn;
    zeroFaceAD = model.AutoDiffBackend.convertToAD(zeros(nf, 1), phi);
    jFaceExternal = zeroFaceAD;
    %jFaceExternal(faces) = -sgn(faces).*current;
    jFaceExternal = subsasgnAD(jFaceExternal,faces, -sgn(faces).*current);
    
    assert(~any(isnan(sgn(faces))));
end
