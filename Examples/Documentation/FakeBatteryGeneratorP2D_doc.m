classdef FakeBatteryGeneratorP2D < BatteryGeneratorP2D


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
            G = computeGeometry(G);

            inputparams.G = G;
            gen.G = G;

        end


    end

    
end

