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

            tbl = crossIndexArray(ytbl, xtbl, {});
            tbl = tbl.addInd('cells', (1 : nx*ny)')';

            celltbl.xind = nxGasSupply + (1 : nxCell)';
            celltbl = IndexArray(celltbl);
            celltbl = crossIndexArray(tbl, celltbl, {'xind'});
            
            params.Cell.cellinds = celltbl.get('cells');

            gassupplytbl.xind = (1 : nxGasSupply)';
            gassupplytbl = IndexArray(gassupplytbl);
            gassupplytbl = crossIndexArray(tbl, gassupplytbl, {'xind'});
                        
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

    end
    
    
end
