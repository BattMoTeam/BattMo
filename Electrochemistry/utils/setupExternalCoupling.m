function [jExternal, jFaceExternal] = setupExternalCoupling(model, phi, phiExternal, conductivity)
    
    coupterm = model.externalCouplingTerm;
    
    jExternal = phi*0.0; %NB hack to initialize zero ad
    
    faces = coupterm.couplingfaces;
    bcval = phiExternal;
    [t, cells] = model.operators.harmFaceBC(conductivity, faces);
    current = t.*(bcval - phi(cells));
    jExternal = subsetPlus(jExternal, current, cells);
    G = model.G;
    nf = G.faces.num;
    sgn = model.operators.sgn;
    zeroFaceAD = model.AutoDiffBackend.convertToAD(zeros(nf, 1), phi);
    jFaceExternal = zeroFaceAD;
    jFaceExternal = subsasgnAD(jFaceExternal, faces, -sgn(faces).*current);
    
    assert(~any(isnan(sgn(faces))));
    
end
