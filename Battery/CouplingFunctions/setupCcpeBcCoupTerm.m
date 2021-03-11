function coupTerm = setupCcpeBcCoupTerm(model)

    ccpe = model.ccpe;
    G = ccpe.G;

    % We pick up the faces at the top of Cccpe
    yf = G.faces.centroids(:, 2);
    myf = max(yf);
    faces = find(abs(yf-myf) < eps*1000);
    cells = sum(G.faces.neighbors(faces, :), 2);

    compnames = {'ccpe'};
    coupTerm = couplingTerm('bc-ccpe', compnames);
    coupTerm.couplingfaces = faces;
    coupTerm.couplingcells = cells;

end
