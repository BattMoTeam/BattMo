classdef ElectrolyserGridGenerator1D < ElectrolyserGridGenerator


    properties

        xs % length of components in x direction
           % x(1) : length of hydrogen porous transport layer
           % x(2) : length of hydrogen catalyst layer
           % x(3) : length of pure ionomer (excluding the catalyst layer)
           % x(4) : length of oxygen catalyst layer
           % x(5) : length of oxygen porous transport layer

        nxs % vector of same length as xs with the corresponding discretization numbers

    end

    methods

        function gen = ElectrolyserGridGenerator1D(xs, nxs)
            gen = gen@ElectrolyserGridGenerator();
            if nargin > 0
                gen.xs = xs;
                gen.nxs = nxs;
            end

        end

        function [paramobj, gen] = updateElectrolyserInputParams(gen, paramobj)

            xs = gen.xs;
            nxs = gen.nxs;

            her = 'HydrogenEvolutionElectrode';
            oer = 'OxygenEvolutionElectrode';
            ptl = 'PorousTransportLayer';
            inm = 'IonomerMembrane';

            params.(inm).cellind = nxs(1) + (1 : sum(nxs(2 : 4)))';
            params.(her).cellind = (1 : sum(nxs(1 : 2)))';
            params.(her).coupcellind = nxs(1) + (1 : nxs(2))'; % coupling with ionomer
            params.(her).bcfaces = 1;
            params.(her).bccells = 1;

            params.(oer).cellind = sum(nxs(1 : 3)) + (1 : sum(nxs(4 : 5)))';
            params.(oer).coupcellind = (1 : nxs(4))';  % coupling with ionomer, note that the indices are given in own grid
            params.(oer).bcfaces = sum(nxs(4 : 5)) + 1; % index in own grid
            params.(oer).bccells = sum(nxs(4 : 5)); % index in own grid

            [paramobj, gen] = gen.setupElectrolyserInputParams(paramobj, params);
        end


        function [paramobj, gen] = setupGrid(gen, paramobj, ~)

            xs = gen.xs;
            nxs = gen.nxs;

            x = xs./nxs;
            x = rldecode(x, nxs);
            x = [0; cumsum(x)];

            G = tensorGrid(x);
            G = computeGeometry(G);

            paramobj.G = G;
            gen.G = G;

        end


    end


end
