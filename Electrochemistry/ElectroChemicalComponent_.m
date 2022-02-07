classdef ElectroChemicalComponent_ < ElectronicComponent_

    methods

        function model = ElectroChemicalComponent_(name)

            model = model@ElectronicComponent_(name);
            
            names = {model.names{:}        , ...
                     'c'
                     'massSource' , ...
                     'massFlux'   , ...
                     'massAccum'  , ...
                     'massCons'};
            model.names = names;
            
            model.vardims('cs') = 2;

            fn = @ElectroChemicalComponent.updateMassFlux;
            inputnames = {'c'};
            fnmodel = {'.'};
            model = model.addPropFunction('massFlux', fn, inputnames, fnmodel);        
            
            fn = @ElectroChemicalComponent.updateMassConservation;
            inputnames = {'massFlux', 'massSource', 'massAccum'};
            fnmodel = {'.'};
            model = model.addPropFunction('massCons', fn, inputnames, fnmodel);

        end
        
    end
    
end




%{
Copyright 2009-2021 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BatMo

BatMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BatMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BatMo.  If not, see <http://www.gnu.org/licenses/>.
%}
