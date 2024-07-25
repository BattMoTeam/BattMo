classdef HomogeneousBlockGridGenerator

    properties

        % Global grid
        parentGrid

        xlength = 1
        ylength = 1
        zlength = 1

        nx = 10
        ny = 10        
        nz = 10

        connectionType = 'lateral faces' % two cases implemented
                                         % - 'lateral faces'
                                         % - 'corners'
        radiusFraction = 0.1 % in case of corner connection
        
    end

    methods

        function [inputparams, gen] = updateGridInputParams(gen, inputparams, params)

            fdnames = {'connectionType'             , ...
                       'radiusFraction'};

            gen = dispatchParams(gen, params, fdnames);
            
            [inputparams, gen] = gen.setupGridInputParams(inputparams, []);

        end


        function [inputparams, gen] = setupGridInputParams(gen, inputparams, params)

            el      = 'ElectronicModel';
            thermal = 'ThermalModel';
            
            [inputparams, gen] = gen.setupGrid(inputparams, params);

            params_el = pickField(params, el);
            inputparams.(el) = gen.setupElectronicModel(inputparams.(el), params_el);

            inputparams = gen.setupExternalElectronicCoupling(inputparams, params);
            
            if inputparams.use_thermal
                
                params_thermal = pickField(params, thermal);
                inputparams.(thermal) = gen.setupThermalModel(inputparams.(thermal), params_thermal);
                inputparams.(thermal) = gen.setupExternalThermalCoupling(inputparams.(thermal), params_thermal);
                
            end

        end


        function [inputparams, gen] = setupGrid(gen, inputparams, params)

            
            G = cartGrid([gen.nx, gen.ny, gen.nz], [gen.xlength, gen.ylength, gen.zlength]);
            
            parentGrid = Grid(G);

            gen.parentGrid = parentGrid;
            
        end


        function inputparams = setupElectronicModel(gen, inputparams, params)

            parentGrid = gen.parentGrid;
            
            inputparams.G = genSubGrid(parentGrid, (1 : parentGrid.getNumberOfCells())');

        end

        function inputparams = setupThermalModel(gen, inputparams, params)

            parentGrid = gen.parentGrid;
            
            inputparams.G = genSubGrid(parentGrid, (1 : parentGrid.getNumberOfCells())');

        end

        function inputparams = setupExternalElectronicCoupling(gen, inputparams, params)

            G = gen.parentGrid.mrstFormat;

            switch gen.connectionType
                
              case 'corners'

                inputcorner  = min(G.nodes.coords, [], 1);
                outputcorner = max(G.nodes.coords, [], 1);

                tbls = setupTables(G, 'includetbls', {'extfacetbl'});
                extfacetbl  = tbls.extfacetbl;
                cellfacetbl = tbls.cellfacetbl;

                if gen.radiusFraction == 0
                    inputcelltbl.cells = 1;
                else
                    d = bsxfun(@minus, G.cells.centroids, inputcorner);
                    d = sum(d.^2, 2);
                    radius = gen.radiusFraction*max([gen.xlength, gen.ylength, gen.zlength]);
                    inputcelltbl.cells = find(d < radius);
                end
                
                inputcelltbl = IndexArray(inputcelltbl);
                
                inputcellfacetbl = crossIndexArray(inputcelltbl, cellfacetbl, {'cells'});
                inputcellfacetbl = crossIndexArray(inputcellfacetbl, extfacetbl, {'faces'});

                coupterm = couplingTerm('input', {'ElectronicModel'});
                coupterm.couplingfaces = inputcellfacetbl.get('faces');
                coupterm.couplingcells = inputcellfacetbl.get('cells');

                couplingTerms{1} = coupterm;

                if gen.radiusFraction == 0
                    outputcelltbl.cells = G.cells.num;
                else
                    d = bsxfun(@minus, G.cells.centroids, outputcorner);
                    d = sum(d.^2, 2);
                    outputcelltbl.cells = find(d < radius);
                end
                
                outputcelltbl = IndexArray(outputcelltbl);
                
                outputcellfacetbl = crossIndexArray(outputcelltbl, cellfacetbl, {'cells'});
                outputcellfacetbl = crossIndexArray(outputcellfacetbl, extfacetbl, {'faces'});

                coupterm = couplingTerm('output', {'ElectronicModel'});
                coupterm.couplingfaces = outputcellfacetbl.get('faces');
                coupterm.couplingcells = outputcellfacetbl.get('cells');

                couplingTerms{2} = coupterm;
                
              case 'lateral faces'

                tbls = setupTables(G);
                cellfacetbl = tbls.cellfacetbl;

                d = max(G.faces.centroids(:, 1)) - min(G.faces.centroids(:, 1));
                tol = 1e-5*d/gen.nx;
                
                inputfacetbl.faces = find(G.faces.centroids(:, 1) < tol);
                inputfacetbl = IndexArray(inputfacetbl);

                inputcellfacetbl = crossIndexArray(inputfacetbl, cellfacetbl, {'faces'});

                coupterm = couplingTerm('input', {'ElectronicModel'});
                coupterm.couplingfaces = inputcellfacetbl.get('faces');
                coupterm.couplingcells = inputcellfacetbl.get('cells');

                couplingTerms{1} = coupterm;

                outputfacetbl.faces = find(G.faces.centroids(:, 1) > max(G.faces.centroids(:, 1)) - tol);
                outputfacetbl = IndexArray(outputfacetbl);

                outputcellfacetbl = crossIndexArray(outputfacetbl, cellfacetbl, {'faces'});

                coupterm = couplingTerm('output', {'ElectronicModel'});
                coupterm.couplingfaces = outputcellfacetbl.get('faces');
                coupterm.couplingcells = outputcellfacetbl.get('cells');

                couplingTerms{2} = coupterm;

              otherwise

                error('connection type not recognized')

            end
            
            inputparams.couplingTerms = couplingTerms;
            
        end

        function inputparams = setupExternalThermalCoupling(gen, inputparams, params)

            G = inputparams.G.mrstFormat;

            tbls = setupTables(G, 'includetbls', {'extfacetbl'});
            extfacetbl  = tbls.extfacetbl;
            cellfacetbl = tbls.cellfacetbl;

            coupcellfacetbl = crossIndexArray(extfacetbl, cellfacetbl, {'faces'});

            coupterm = couplingTerm('external', {'External'});
            coupterm.couplingfaces = coupcellfacetbl.get('faces');
            coupterm.couplingcells = coupcellfacetbl.get('cells');

            inputparams.couplingTerm = coupterm;
            
        end
        
    end


end


