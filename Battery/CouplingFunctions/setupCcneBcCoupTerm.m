function coupTerm = setupCcneBcCoupTerm(model)

    ccne = model.ccne;
    G = ccne.G;

    % We pick up the faces at the top of Cccne
    yf = G.faces.centroids(:, 2);
    myf = max(yf);
    faces = find(abs(yf-myf) < eps*1000);
    cells = sum(G.faces.neighbors(faces, :), 2);

    compnames = {'ccne'};
    coupTerm = couplingTerm('bc-ccne', compnames);
    coupTerm.couplingfaces = faces;
    coupTerm.couplingcells = cells;

end
