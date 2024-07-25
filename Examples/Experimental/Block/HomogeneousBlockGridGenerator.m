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
        
    end

    methods

        function [inputparams, gen] = updateGridInputParams(gen, inputparams)

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
            
            tbls = setupTables(G);
            cellfacetbl = tbls.cellfacetbl;
            
            inputfacetbl.faces = find(G.faces.centroids(:, 1) == 0);
            inputfacetbl = IndexArray(inputfacetbl);

            inputcellfacetbl = crossIndexArray(inputfacetbl, cellfacetbl, {'faces'});

            coupterm = couplingTerm('input', {'ElectronicModel'});
            coupterm.couplingfaces = inputcellfacetbl.get('faces');
            coupterm.couplingcells = inputcellfacetbl.get('cells');

            couplingTerms{1} = coupterm;

            outputfacetbl.faces = find(G.faces.centroids(:, 1) == max(G.faces.centroids(:, 1)));
            outputfacetbl = IndexArray(outputfacetbl);

            outputcellfacetbl = crossIndexArray(outputfacetbl, cellfacetbl, {'faces'});

            coupterm = couplingTerm('output', {'ElectronicModel'});
            coupterm.couplingfaces = outputcellfacetbl.get('faces');
            coupterm.couplingcells = outputcellfacetbl.get('cells');

            couplingTerms{2} = coupterm;

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


