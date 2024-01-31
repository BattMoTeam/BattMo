classdef SectorBatteryGenerator < SpiralBatteryGenerator

    methods

        function gen = SectorBatteryGenerator()
            gen = gen@SpiralBatteryGenerator();
        end

        function [inputparams, gen] = updateBatteryInputParams(gen, inputparams, params)

            gen.nwindings    = params.nwindings;
            gen.rInner       = params.rInner;
            gen.widthDict    = params.widthDict;
            gen.nrDict       = params.nrDict;
            gen.nas          = params.nas;
            gen.L            = params.L;
            gen.nL           = params.nL;

            gen.use_thermal = inputparams.use_thermal;

            [inputparams, gen] = gen.setupBatteryInputParams(inputparams, []);

        end

        function [inputparams, gen] = setupGrid(gen, inputparams, ~)

            [gen, G] = sectorGrid(gen);
            inputparams.G = G;

        end

    end

end


%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
