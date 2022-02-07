classdef Electrolyte_ < ElectroChemicalComponent_
%
% The function that couples the electolyte and the electrode is not implemented here (but at the battery level)
%
%
    methods
        
        function model = Electrolyte_(name)
            
            model = model@ElectroChemicalComponent_(name);
            
            names = {model.names{:} , ...
                     'cs'            , ...
                     'D'            , ...
                     'dmudcs'       , ...
                     'conductivity' , ...
                     'jchems'       , ...
                     'diffFlux'};
            model.names = names;
            
            model.vardims('cs') = 2;
            model.vardims('dmudcs') = 2;
            model.vardims('jchems') = 2;
            
            fn = @Electrolyte.updateConductivity;
            inputnames = {'c'}
            fnmodel = {'.'};
            model = model.addPropFunction('cs', fn, inputnames, fnmodel);
            
            fn = @Electrolyte.updateConductivity;
            inputnames = {'c', 'T', 'phi'};
            fnmodel = {'.'};
            model = model.addPropFunction('conductivity', fn, inputnames, fnmodel);

            
            fn = @Electrolyte.updateChemicalCurrent;
            inputnames = {'c', 'T', 'phi'};
            fnmodel = {'.'};
            model = model.addPropFunction('dmudcs', fn, inputnames, fnmodel);
            model = model.addPropFunction('jchems', fn, inputnames, fnmodel);
            
            fn = @Electrolyte.updateDiffusionCoefficient;
            inputnames = {'c', 'T'};
            fnmodel = {'.'};
            model = model.addPropFunction('D', fn, inputnames, fnmodel);

            fn = @Electrolyte.updateCurrent;
            inputnames = {'phi', 'jchems', 'conductivity'};
            fnmodel = {'.'};
            model = model.addPropFunction('j', fn, inputnames, fnmodel);
            
            fn = @Electrolyte.updateMassFlux;
            inputnames = {'c', 'j', 'D'};
            fnmodel = {'.'};
            model = model.addPropFunction('massFlux', fn, inputnames, fnmodel);
            model = model.addPropFunction('diffFlux', fn, inputnames, fnmodel);

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
