classdef FakeBatteryGeneratorP2D_doc < BatteryGeneratorP2D


    methods

        function [inputparams, gen] = setupGrid(gen, inputparams, ~)

            sepnx  = gen.sepnx;
            nenx   = gen.nenx;
            penx   = gen.penx;
            ccnenx = gen.ccnenx;
            ccpenx = gen.ccpenx;

            xlength = gen.xlength;

            if gen.include_current_collectors
                nxs = [ccnenx; nenx; sepnx; penx; ccpenx];
            else
                nxs = [nenx; sepnx; penx];
                xlength = xlength(2 : 4);
            end

            x = xlength./nxs;
            x = rldecode(x, nxs);
            x = [0; cumsum(x)];

            d = sqrt(gen.faceArea);
            y = [0; d];
            z = [0; d];
            
            G = tensorGrid(x, y, z);

            parentGrid = Grid(G, 'faceArea', gen.faceArea);

            G = genSubGrid(parentGrid, (1 : parentGrid.getNumberOfCells())');

            inputparams.G  = G;
            gen.parentGrid = parentGrid;


        end


    end

    
end

