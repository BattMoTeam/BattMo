function [output, G] = sectorGrid(params)

    %% Parameters
    %
    % params structure with following fields
    % - nwindings : number of windings in the spiral
    % - rInner        : "radius" at the middle
    % - widthDict : dictionary of widths for each component. The required key names for the dictionary are
    %                 - 'Separator'
    %                 - 'NegativeActiveMaterial'
    %                 - 'NegativeCurrentCollector'
    %                 - 'PositiveActiveMaterial'
    %                 - 'PositiveCurrentCollector'
    % - nrDict    : dicionary with number of cell in radial direction for each component (same keys as in widthDict).
    % - L         : length of the battery
    % - nas       : number of cells in the angular direction
    % - nL        : number of discretization cells in the longitudonal
    %
    % RETURNS
    %
    % - G       : grid
    % - tag     : cell-valued vector giving component number (indexing is given by tagdict)
    % - tagdict : dictionary giving the component number


    output = radialGrid(params);

    nwindings = params.nwindings;
    nrDict    = params.nrDict;
    nas       = params.nas;
    nL        = params.nL;

    nrs = [nrDict('PositiveCoating')          ; ...
           nrDict('PositiveCurrentCollector') ; ...
           nrDict('PositiveCoating')          ; ...
           nrDict('Separator')                ; ...
           nrDict('NegativeCoating')          ; ...
           nrDict('NegativeCurrentCollector') ; ...
           nrDict('NegativeCoating')          ; ...
           nrDict('Separator')];

    nR = sum(nrs)*nwindings;

    radG                    = output.parentGrid;
    tag                     = output.tag;
    tagdict                 = output.tagdict;
    positiveExtCurrentFaces = output.positiveExtCurrentFaces;
    negativeExtCurrentFaces = output.negativeExtCurrentFaces;
    thermalExchangeFaces    = output.thermalExchangeFaces;

    [indA, indR, indZ] = ind2sub([nas, nR, nL], (1 : radG.getNumberOfCells())');
    celltbl.cells = (1 : radG.getNumberOfCells())';
    celltbl.indA = indA;
    celltbl.indR = indR;
    celltbl.indZ = indZ;
    celltbl = IndexArray(celltbl);

    scelltbl.indA = 1;
    scelltbl = IndexArray(scelltbl);

    scelltbl = crossIndexArray(scelltbl, celltbl, {'indA'});

    scells = scelltbl.get('cells');

    rcells = (1 : radG.getNumberOfCells())';
    rcells(scells) = [];

    mradG = radG.mrstFormat();
    [G, cellmap, facemap] = removeCells(mradG, rcells);

    cutcelltbl.radcells = cellmap;
    cutcelltbl.cells    = (1 : G.cells.num)';
    cutcelltbl = IndexArray(cutcelltbl);

    gen = CrossIndexArrayGenerator();
    gen.tbl1 = cutcelltbl;
    gen.tbl2 = celltbl;
    gen.replacefds2 = {{'cells', 'radcells'}};
    gen.mergefds = {'radcells'};

    celltbl = gen.eval();
    celltbl = projIndexArray(celltbl, {'cells', 'indR', 'indZ'});

    tag = tag(cellmap);

    tbls = setupTables(mradG);
    radfacetbl = tbls.facetbl;
    tbls = setupTables(G);
    facetbl = tbls.facetbl;
    facetbl = facetbl.addInd('radfaces', facemap);

    map = TensorMap();
    map.fromTbl = radfacetbl;
    map.toTbl = facetbl;
    map.replaceFromTblfds = {{'faces', 'radfaces'}};
    map.mergefds = {'radfaces'};
    map = map.setup();

    sfaces = zeros(radfacetbl.num, 1);
    sfaces(positiveExtCurrentFaces) = 1;
    sfaces = map.eval(sfaces);
    positiveExtCurrentFaces = find(sfaces);

    sfaces = zeros(radfacetbl.num, 1);
    sfaces(negativeExtCurrentFaces) = 1;
    sfaces = map.eval(sfaces);
    negativeExtCurrentFaces = find(sfaces);

    sfaces = zeros(radfacetbl.num, 1);
    sfaces(thermalExchangeFaces) = 1;
    sfaces = map.eval(sfaces);
    thermalExchangeFaces = find(sfaces);

    parentGrid = Grid(G);

    G = genSubGrid(parentGrid, (1 : parentGrid.getNumberOfCells())');

    %% setup output structure
    output = params;
    output.parentGrid              = parentGrid;
    output.tag                     = tag;
    output.tagdict                 = tagdict;
    output.celltbl                 = celltbl;
    output.positiveExtCurrentFaces = positiveExtCurrentFaces;
    output.negativeExtCurrentFaces = negativeExtCurrentFaces;
    output.thermalExchangeFaces    = thermalExchangeFaces;

end


%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
