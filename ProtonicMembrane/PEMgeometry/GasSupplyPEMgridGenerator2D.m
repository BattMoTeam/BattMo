classdef GasSupplyPEMgridGenerator2D < GasSupplyPEMgridGenerator

    properties

        nxCell
        nxGasSupply
        
        lxCell
        lxGasSupply

        nx
        ny
        lx
        ly

    end

    methods
        
        function [paramobj, gen] = updateInputParams(gen, paramobj, params)

            nxCell      = gen.nxCell;
            nxGasSupply = gen.nxGasSupply;
            ny          = gen.ny;
            
            nx = nxCell + nxGasSupply;
            lx = gen.lxCell + gen.lxGasSupply;

            gen.nx = nx;
            gen.lx = lx;
            
            xtbl.xind = (1 : nx)';
            xtbl = IndexArray(xtbl);
            ytbl.yind = (1 : ny)';
            ytbl = IndexArray(ytbl);

            carttbl = crossIndexArray(ytbl, xtbl, {});
            carttbl = carttbl.addInd('cells', (1 : nx*ny)')';

            celltbl.xind = nxGasSupply + (1 : nxCell)';
            celltbl = IndexArray(celltbl);
            celltbl = crossIndexArray(carttbl, celltbl, {'xind'});
            
            params.Cell.cellinds = celltbl.get('cells');

            gassupplytbl.xind = (1 : nxGasSupply)';
            gassupplytbl = IndexArray(gassupplytbl);
            gassupplytbl = crossIndexArray(carttbl, gassupplytbl, {'xind'});
                        
            params.GasSupply.cellinds = gassupplytbl.get('cells');
            
            [paramobj, gen] = setupInputParams(gen, paramobj, params);
            
        end

        function [paramobj, gen] = setupGrid(gen, paramobj, params)
            
            nx1 = gen.nxGasSupply;
            nx2 = gen.nxCell;
            lx1 = gen.lxGasSupply;
            lx2 = gen.lxCell;
            ny  = gen.ny;
            ly  = gen.ly;

            x = [(lx1/nx1)*ones(nx1, 1); (lx2/nx2)*ones(nx2, 1)];
            x = [0; cumsum(x)];
            y = (ly/ny)*ones(ny, 1);
            y = [0; cumsum(y)];

            G = tensorGrid(x, y);

            paramobj.G = G;
            gen.G      = G;
            
        end
        
        function [paramobj, gen] = setupCellGridGenerator(gen, paramobj, params)
        % setup cellGridGenerator

            cgen = PEMgridGenerator2D();

            cgen.xlength = gen.lxCell;
            cgen.ylength = gen.ly;
            cgen.Nx      = gen.nxCell;
            cgen.Ny      = gen.ny;

            gen.cellGridGenerator = cgen;
            
        end

        function [paramobj, gen] = setupGasSupplyGridGenerator(gen, paramobj, params)
        % setup gasSupplyGridGenerator

            gsgen = GasSupplyGridGenerator2D();

            gsgen.nx = gen.nxGasSupply;
            gsgen.ny = gen.ny;
            gsgen.lx = gen.lxGasSupply;
            gsgen.ly = gen.ly;

            gen.gasSupplyGridGenerator = gsgen;
            
        end

        function [paramobj, gen] = setupCellGasSupplyCoupling(gen, paramobj, params)


            coupTerms = paramobj.Cell.couplingTerms;

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

            % We setup Cell mapping
            cG = paramobj.Cell.G;
            cfacefacetbl.cfaces = (1 : cG.faces.num)';
            cfacefacetbl.faces = cG.mappings.facemap;
            cfacefacetbl = IndexArray(cfacefacetbl);
            % We setup GasSupply mapping
            gsG = paramobj.GasSupply.G;
            gfacefacetbl.gfaces = (1 : gsG.faces.num)';
            gfacefacetbl.faces = gsG.mappings.facemap;
            gfacefacetbl = IndexArray(gfacefacetbl);

            coupcfacefacetbl      = crossIndexArray(coupcfacetbl, cfacefacetbl, {'cfaces'});
            coupgfacecfacefacetbl = crossIndexArray(coupcfacefacetbl, gfacefacetbl, {'faces'});

            tbls = setupSimpleTables(gsG);
            gcellgfacetbl = tbls.cellfacetbl;
            gcellgfacetbl = replacefield(gcellgfacetbl, {{'cells', 'gcells'}, {'faces', 'gfaces'}});
            

            coupgcellgfacecfacefacetbl = crossIndexArray(coupgfacecfacefacetbl, gcellgfacetbl, {'gfaces'});

            coupTerm = couplingTerm('Gas Supply - Anode', {'GasSupply', {'Cell', 'Anode'}});
            tbl = coupgcellgfacecfacefacetbl; % shortcut
            coupTerm.couplingcells = [tbl.get('gcells'), tbl.get('ind')];
            coupTerm.couplingfaces = [tbl.get('gfaces'), tbl.get('ind')];

            paramobj.couplingTerm = coupTerm;
            
        end        
    end
    
    
end
