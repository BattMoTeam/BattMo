classdef CellGridGenerator2D < CellGridGenerator

    properties

        nxElectrolyser
        nxGasSupply
        
        lxElectrolyser
        lxGasSupply

        nx
        ny
        lx
        ly

    end

    methods
        
        function [inputparams, gen] = updateInputParams(gen, inputparams, params)

            % We start with grid setup
            nx1 = gen.nxGasSupply;
            nx2 = gen.nxElectrolyser;
            lx1 = gen.lxGasSupply;
            lx2 = gen.lxElectrolyser;
            ny  = gen.ny;
            ly  = gen.ly;

            x = [(lx1/nx1)*ones(nx1, 1); (lx2/nx2)*ones(nx2, 1)];
            x = [0; cumsum(x)];
            y = (ly/ny)*ones(ny, 1);
            y = [0; cumsum(y)];

            G = tensorGrid(x, y);

            parentGrid = Grid(G);

            G = genSubGrid(parentGrid, (1 : parentGrid.getNumberOfCells())');
            
            inputparams.G = G;
            gen.parentGrid = parentGrid;

            % We setup params structure which defines which cells are in the PEM and GasSupply
            % NOTE that we rely implicity below on ordering of cells used in tensorGrid
            
            nxElectrolyser = gen.nxElectrolyser;
            nxGasSupply    = gen.nxGasSupply;
            ny             = gen.ny;
            
            nx = nxElectrolyser + nxGasSupply;
            lx = gen.lxElectrolyser + gen.lxGasSupply;

            gen.nx = nx;
            gen.lx = lx;
            
            xtbl.xind = (1 : nx)';
            xtbl = IndexArray(xtbl);
            ytbl.yind = (1 : ny)';
            ytbl = IndexArray(ytbl);

            carttbl = crossIndexArray(ytbl, xtbl, {});
            % NOTE that we rely implicity on ordering of cells used in tensorGrid
            carttbl = carttbl.addInd('cells', (1 : nx*ny)')';

            celltbl.xind = nxGasSupply + (1 : nxElectrolyser)';
            celltbl = IndexArray(celltbl);
            celltbl = crossIndexArray(carttbl, celltbl, {'xind'});
            
            params.Electrolyser.cellinds             = celltbl.get('cells');
            params.Electrolyser.Electrolyte.cellinds = params.Electrolyser.cellinds;
            
            gassupplytbl.xind = (1 : nxGasSupply)';
            gassupplytbl = IndexArray(gassupplytbl);
            gassupplytbl = crossIndexArray(carttbl, gassupplytbl, {'xind'});
                        
            params.GasSupply.cellinds = gassupplytbl.get('cells');
            
            [inputparams, gen] = setupInputParams(gen, inputparams, params);
            
        end

        function [inputparams, gen] = setupGrid(gen, inputparams, params)

        % nothing is done here. The grid is setup in the updateInputParams method.
            
        end
        
        function gen = setupElectrolyserGridGenerator(gen, inputparams, params)
        % setup electrolyserGridGenerator

            cgen = PEMgridGenerator2D();

            cgen.parentGrid = gen.parentGrid;
            cgen.xlength    = gen.lxElectrolyser;
            cgen.ylength    = gen.ly;
            cgen.Nx         = gen.nxElectrolyser;
            cgen.Ny         = gen.ny;
            
            gen.electrolyserGridGenerator = cgen;
            
        end

        function gen = setupGasSupplyGridGenerator(gen, inputparams, params)
        % setup gasSupplyGridGenerator

            gsgen = GasSupplyGridGenerator2D();

            gsgen.parentGrid = gen.parentGrid;
            gsgen.Nx         = gen.nxGasSupply;
            gsgen.Ny         = gen.ny;
            gsgen.xlength    = gen.lxGasSupply;
            gsgen.ylength    = gen.ly;

            gen.gasSupplyGridGenerator = gsgen;
            
        end

        function [inputparams, gen] = setupElectrolyserGasSupplyCoupling(gen, inputparams, params)


            coupTerms = inputparams.Electrolyser.couplingTerms;

            for icoup = 1 : numel(coupTerms)
                coupTerm = coupTerms{icoup};
                [lia, locb] = ismember('Anode', coupTerm.componentnames);
                if lia
                    loc = true(2, 1);
                    loc(locb) = false;
                    cfaces = coupTerm.couplingfaces(:, loc);
                    break
                end
            end

            coupcfacetbl.cfaces = cfaces;
            coupcfacetbl.ind = (1 : numel(cfaces))'; % we keep track of local indexing which is used for the anode structure
            coupcfacetbl = IndexArray(coupcfacetbl);

            % We setup Electrolyser mapping
            cG = inputparams.Electrolyser.G;
            cfacefacetbl.cfaces = (1 : cG.getNumberOfFaces())';
            cfacefacetbl.faces = cG.mappings.facemap;
            cfacefacetbl = IndexArray(cfacefacetbl);
            % We setup GasSupply mapping
            gsG = inputparams.GasSupply.G;
            gfacefacetbl.gfaces = (1 : gsG.getNumberOfFaces())';
            gfacefacetbl.faces = gsG.mappings.facemap;
            gfacefacetbl = IndexArray(gfacefacetbl);

            coupcfacefacetbl      = crossIndexArray(coupcfacetbl, cfacefacetbl, {'cfaces'});
            coupgfacecfacefacetbl = crossIndexArray(coupcfacefacetbl, gfacefacetbl, {'faces'});

            tbls = setupTables(gsG.mrstFormat);
            gcellgfacetbl = tbls.cellfacetbl;
            gcellgfacetbl = replacefield(gcellgfacetbl, {{'cells', 'gcells'}, {'faces', 'gfaces'}});
            

            coupgcellgfacecfacefacetbl = crossIndexArray(coupgfacecfacefacetbl, gcellgfacetbl, {'gfaces'});

            coupTerm = couplingTerm('Gas Supply - Anode', {'GasSupply', {'Electrolyser', 'Anode'}});
            tbl = coupgcellgfacecfacefacetbl; % shortcut
            coupTerm.couplingcells = [tbl.get('gcells'), tbl.get('ind')];
            coupTerm.couplingfaces = [tbl.get('gfaces'), tbl.get('ind')];

            inputparams.couplingTerm = coupTerm;
            
        end        
    end
    
    
end
