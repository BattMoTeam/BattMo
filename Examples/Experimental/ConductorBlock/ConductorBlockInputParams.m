classdef ConductorBlockInputParams < InputParams
%
% Input parameter class for the :code:`Battery` model.
%

    properties

        ElectronicModel
        ThermalModel        
        Control
        
        couplingTerms % Coupling terms
        use_thermal

        initT % Initial temperature

    end

    methods

        function inputparams = ConductorBlockInputParams(jsonstruct)

            el      = 'ElectronicModel';
            thermal = 'ThermalModel';
            ctrl    = 'Control';
            
            jsonstruct = equalizeJsonStructField(jsonstruct, 'use_thermal', {el, 'use_thermal'});

            inputparams = inputparams@InputParams(jsonstruct);

            pick = @(fd) pickField(jsonstruct, fd);

            inputparams.(el)      = ElectronicComponentInputParams(pick(el));
            inputparams.(thermal) = ThermalComponentInputParams(pick(thermal));
            inputparams.(ctrl)    = ConductorBlockControlModelInputParams(pick(ctrl));

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
