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

            nxCell     = gen.nxCell;
            nxGasSupply = gen.nxGasSupply;
            ny         = gen.ny;
            
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

            gaslayertbl.xind = (1 : nxGasSupply)';
            gaslayertbl = IndexArray(gaslayertbl);
            gaslayertbl = crossIndexArray(tbl, gaslayertbl, {'xind'});
                        
            params.GasSupply.cellinds = gaslayertbl.get('cells');
            
            [paramobj, gen] = setupInputParams(gen, paramobj, params);
            
        end

        function [paramobj, gen] = setupGrid(gen, paramobj, params)
            
            nx = gen.nx;
            ny = gen.ny;
            lx = gen.lx;
            ly = gen.ly;

            G = cartGrid([nx, ny], [lx, ly]);

            paramobj.G = G;
            gen.G      = G;
            
        end
        
        function [paramobj, gen] = setupCellGridGenerator(gen, paramobj, params)
        % setup cellGridGenerator

            gen.cellGridGenerator = PEMgridGenerator2D();
            
        end

        function [paramobj, gen] = setupGasSupplyGridGenerator(gen, paramobj, params)
        % setup gasLayerGridGenerator
            
            gen.gasLayerGridGenerator = GasSupplyGridGenerator2D();
            
        end

    end
    
    
end
