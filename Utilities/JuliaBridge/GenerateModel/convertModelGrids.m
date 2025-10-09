function model = convertModelGrids(model)
% Convert the grids from sub-grid representation to a concrete grid with the data needed for a julia simulation.
    
    if isprop(model, 'G') && ~isempty(model.G)

        %% We need to provide the following
        % volumes        % volumes
	    % neighborship   % Internal faces only
        % half_trans     % half transmissibilities for the internal faces
	    % boundary_cells % indices of the boundary cells (some can can be repeated if a cell has two boundary faces). Same length as boundary_hfT.
	    % boundary_hfT   % Boundary half face transmissibilities

        G = model.G;
        
        volumes =  G.getVolumes();
        
        tbls = setupTables(G.mrstFormat, 'includetbls', {'intfacetbl', 'cellintfacetbl', 'extfacetbl'});

        cellintfacetbl = tbls.cellintfacetbl;
        cellfacetbl    = tbls.cellfacetbl;
        intfacetbl     = tbls.intfacetbl;
        extfacetbl     = tbls.extfacetbl;
        
        tpfvGeometry = G.getTPFVgeometry();
        hT           = tpfvGeometry.hT; % hT is in cellfacetbl

        %% Setup internal half transmissibilities

        cellintfacetbl = sortIndexArray(cellintfacetbl, {'faces', 'cells'});
        neighborship = cellintfacetbl.get('cells');
        neighborship = reshape(neighborship, 2, []);

        map = TensorMap();
        map.fromTbl  = cellfacetbl;
        map.toTbl    = cellintfacetbl;
        map.mergefds = {'cells', 'faces'};
        map = map.setup();

        half_trans = map.eval(hT);
        half_trans = reshape(half_trans, 2, []);
        
        %% Setup external half transmissibilities
        
        cellextfacetbl = crossIndexArray(cellfacetbl, extfacetbl, {'faces'});

        boundary_cells = cellextfacetbl.get('cells');

        map = TensorMap();
        map.fromTbl  = cellfacetbl;
        map.toTbl    = cellextfacetbl;
        map.mergefds = {'cells', 'faces'};
        map = map.setup();

        boundary_hfT = map.eval(hT);
        
        model.G = struct('volumes'       , volumes       , ...
	                     'neighborship'  , neighborship  , ...
                         'half_trans'    , half_trans    , ...
	                     'boundary_cells', boundary_cells, ...
	                     'boundary_hfT'  , boundary_hfT);

    end

    submodelnames = model.getSubModelNames();

    if ~isempty(submodelnames)

        for isub = 1 : numel(submodelnames)

            submodelname = submodelnames{isub};

            model.(submodelname) = convertModelGrids(model.(submodelname));

        end

    end
    
end

