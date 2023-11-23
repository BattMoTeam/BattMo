classdef GasSupplyGridGenerator2D < GasSupplyGridGenerator
    
    properties
        
        nx
        ny
        lx
        ly

    end


    methods

        
        function [paramobj, gen] = updateInputParams(gen, paramobj, params)

            paramobj = gen.setupInputParams(paramobj, []);
            
        end
        

        function [paramobj, gen] = setupGrid(gen, paramobj, params)

            nx = gen.nx;
            ny = gen.ny;
            lx = gen.lx;
            ly = gen.ly;

            G = cartGrid([nx, ny], [lx, ly]);
            G = computeGeometry(G);

            paramobj.G = G;
            gen.G      = G;
            
        end


        function [paramobj, gen] = setupExternalCoupling(gen, paramobj, params)

            G  = gen.G;
            nx = gen.nx;
            ny = gen.ny;
            lx = gen.lx;
            ly = gen.ly;
            
            tbls = setupSimpleTables(G);
            cellfacetbl = tbls.cellfacetbl;

            clear bcfacecouptbl1;
            bcfacecouptbl1.faces = (nx + 1)*ny + (1 : floor(nx/2))';
            bcfacecouptbl1 = IndexArray(bcfacecouptbl1);
            nc = bcfacecouptbl1.num;
            bcfacecouptbl1 =  bcfacecouptbl1.addInd('coup', ones(nc, 1));

            clear bcfacecouptbl2;
            bcfacecouptbl2.faces = (nx + 1)*ny + ny*nx + (floor(nx/2) + 1 : nx)';
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

            paramobj.couplingTerms = couplingTerms;

        end

    end

end
