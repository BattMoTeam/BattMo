classdef LengthSetter1D

    properties

        nxs        % number of cells for each component
        reflengths % reference lengths for each component
        compinds   % index of the components

    end

    methods

        function lengthsetter = LengthSetter1D(gridGenerator, components)
        %  gridGenerator is the BatteryGenerator1D structure
            assert(~gridGenerator.include_current_collectors, 'we do not include current collectors')

            gen = gridGenerator; % shortcut

            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            sep = 'Separator';

            compinds = ismember({ne, sep, pe}, components);

            lengthsetter.reflengths = gen.xlength(2 : 4);
            lengthsetter.nxs        = [gen.nenx; gen.sepnx; gen.penx];
            lengthsetter.compinds   = compinds;

        end

        function model = setLengths(lengthsetter, model, v)

            nxs        = lengthsetter.nxs;
            reflengths = lengthsetter.reflengths;
            compinds   = lengthsetter.compinds;

            lgths = reflengths;
            lgths = subsasgnAD(lgths, compinds, v);

            nx = sum(nxs);

            x = zeros(nx + 1, 1);

            x(1) = 0;
            startind = 1;
            for i = 1 : 3
                dx = lgths(i)/nxs(i).*ones(nxs(i), 1);
                x = subsasgnAD(x, startind + (1 : nxs(i)), x(startind) + cumsum(dx));
                startind = startind + nxs(i);
            end

            [tPFVgeometry, model.G.parentGrid] = model.G.parentGrid.computeTPFVgeometry(x, 'assemblyType', 'AD');
            model = model.setTPFVgeometry(tPFVgeometry);

        end

        function v = getAllLengths(lengthsetter, model)
        % v is length value and comp is one component name

            nxs      = lengthsetter.nxs;

            x = model.G.parentGrid.tPFVgeometry.nodes.coords;

            nxs = cumsum(nxs);

            v = double2ADI(zeros(2, 1), x);
            xstart = 0;
            for i = 1 : 3
                v(i) = x(nxs(i) + 1) - xstart;
                xstart = x(nxs(i) + 1);
            end

        end


        function v = getLengths(lengthsetter, model)

            compinds = lengthsetter.compinds;

            v = lengthsetter.getAllLengths(model);
            v = v(compinds);

        end
    end

end
