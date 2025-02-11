classdef GasSupplyGridGenerator2D < GasSupplyGridGenerator
    
    properties
        
        Nx
        Ny
        xlength
        ylength

    end


    methods

        
        function [inputparams, gen] = updateInputParams(gen, inputparams, params)

            inputparams = gen.setupInputParams(inputparams, []);
            
        end
        

        function [inputparams, gen] = setupGrid(gen, inputparams, params)

            Nx      = gen.Nx;
            Ny      = gen.Ny;
            xlength = gen.xlength;
            ylength = gen.ylength;

            G = cartGrid([Nx, Ny], [xlength, ylength]);
            
            parentGrid = Grid(G);
            
            G = genSubGrid(parentGrid, (1 : parentGrid.getNumberOfCells())');
            
            inputparams.G = G;
            gen.parentGrid = parentGrid;
            
        end


        function inputparams = setupExternalCoupling(gen, inputparams, params)

            G  = inputparams.G.mrstFormat;
            
            Nx      = gen.Nx;
            Ny      = gen.Ny;
            xlength = gen.xlength;
            ylength = gen.ylength;
            
            tbls = setupTables(G);
            cellfacetbl = tbls.cellfacetbl;

            clear bcfacecouptbl1;
            bcfacecouptbl1.faces = (Nx + 1)*Ny + (1 : floor(Nx/2))';
            bcfacecouptbl1 = IndexArray(bcfacecouptbl1);
            nc = bcfacecouptbl1.num;
            bcfacecouptbl1 =  bcfacecouptbl1.addInd('coup', ones(nc, 1));

            clear bcfacecouptbl2;
            bcfacecouptbl2.faces = (Nx + 1)*Ny + Ny*Nx + (1 : floor(Nx/2))';
            bcfacecouptbl2 = IndexArray(bcfacecouptbl2);
            nc = bcfacecouptbl2.num;
            bcfacecouptbl2 =  bcfacecouptbl2.addInd('coup', 2*ones(nc, 1));

            bcfacecouptbl = concatIndexArray(bcfacecouptbl1, bcfacecouptbl2, {});

            bccellfacecouptbl = crossIndexArray(bcfacecouptbl, cellfacetbl, {'faces'});

            coupnames = {'input', 'output'};

            couplingTerms = {};

            for  icoup = 1 : numel(coupnames)

                coupname = coupnames{icoup};
                name = sprintf('External %s', coupname);
                compnames = {'External'};
                coupTerm = couplingTerm(name, compnames);

                clear couptbl;
                couptbl.coup =  icoup;
                couptbl = IndexArray(couptbl);

                tbl =  crossIndexArray(couptbl, bccellfacecouptbl, {'coup'});

                coupTerm.couplingfaces = tbl.get('faces');
                coupTerm.couplingcells = tbl.get('cells');

                couplingTerms{end + 1} = coupTerm;
                
            end

            inputparams.couplingTerms = couplingTerms;

        end

    end

end
