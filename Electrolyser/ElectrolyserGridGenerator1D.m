classdef ElectrolyserGridGenerator1D < ElectrolyserGridGenerator


    properties

        xs % length of components in x direction
           % x(1) : length of oxygen porous transport layer
           % x(2) : length of oxygen catalyst layer
           % x(3) : length of pure ionomer (excluding the catalyst layer)
           % x(4) : length of hydrogen catalyst layer
           % x(5) : length of hydrogen porous transport layer

        cxs % length of the finite volume cells for each component (same dimension as )
        nxs % vector of same length as xs with the corresponding discretization numbers

    end

    methods

        function gen = ElectrolyserGridGenerator1D(xs, cxs)
            gen = gen@ElectrolyserGridGenerator();
            if nargin > 0
                gen.xs = xs;
                gen.cxs = cxs;
            end

        end

        function [paramobj, gen] = updateElectrolyserInputParams(gen, paramobj)

            xs  = gen.xs;
            cxs = gen.cxs;

            nxs = round(xs./cxs);
            gen.nxs = nxs;

            oer = 'OxygenEvolutionElectrode';
            her = 'HydrogenEvolutionElectrode';
            ptl = 'PorousTransportLayer';
            ctl = 'CatalystLayer';
            inm = 'IonomerMembrane';

            params.(inm).cellind = nxs(1) + (1 : sum(nxs(2 : 4)))';

            params.(oer).cellind       = (1 : sum(nxs(1 : 2)))';
            params.(oer).(ptl).cellind = (1 : sum(nxs(1 : 2)))';
            params.(oer).(ctl).cellind = nxs(1) + (1 : nxs(2))';
            params.(oer).bcfaces = 1;
            params.(oer).bccells = 1;

            params.(her).cellind       = sum(nxs(1 : 3)) + (1 : sum(nxs(4 : 5)))';
            params.(her).(ptl).cellind = sum(nxs(1 : 3)) + (1 : sum(nxs(4 : 5)))';
            params.(her).(ctl).cellind = sum(nxs(1 : 3)) + (1 : nxs(4))';
            params.(her).bcfaces = sum(nxs(4 : 5)) + 1; % index in own grid
            params.(her).bccells = sum(nxs(4 : 5)); % index in own grid

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
