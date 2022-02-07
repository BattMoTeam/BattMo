classdef ElectronicComponent_ < CompositeModel

    methods

        function model = ElectronicComponent_(name)

            model = model@CompositeModel(name);
            names = {'T'        , ...
                     'phi'      , ...
                     'jBcSource', ...
                     'eSource'  , ...
                     'j'        , ...
                     'chargeCons'};
            model.names = names;
            
            fn = @ElectronicComponent.updateCurrent;
            inputnames = {'phi'};
            fnmodel = {'.'};
            model = model.addPropFunction('j', fn, inputnames, fnmodel);
            
            fn = @ElectronicComponent.updateChargeConservation;
            inputnames = {'j', 'jBcSource', 'eSource'};
            fnmodel = {'.'};
            model = model.addPropFunction('chargeCons', fn, inputnames, fnmodel);
            
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
