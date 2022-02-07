classdef ElectrodeActiveComponent_ < ElectroChemicalComponent_

    methods

        function model = ElectrodeActiveComponent_(name)

            model = model@ElectroChemicalComponent_(name);
            
            names = {model.names{:}, ...
                     'jCoupling'};
            model.names = names;

            model.SubModels{1} = ActiveMaterial_('am');
            
            fn = @ElectrodeActiveComponent.updatejBcSource;
            inputnames = {'jCoupling'};
            fnmodel = {'.'};
            model = model.addPropFunction('jBcSource', fn, inputnames, fnmodel);            

            fn = @ElectrodeActiveComponent.updateIonAndCurrentSource;
            inputnames = {VarName({'am'}, 'R')};
            fnmodel = {'.'};
            model = model.addPropFunction('massSource', fn, inputnames, fnmodel);
            model = model.addPropFunction('eSource', fn, inputnames, fnmodel);

            fn = @ElectrodeActiveComponent.updateChargeCarrier;
            inputnames = {VarName({'..'}, 'c')};
            fnmodel = {'..'};
            model = model.addPropFunction({'am', 'cElectrodeAveraged'}, fn, inputnames, fnmodel);

            fn = @ElectrodeActiveComponent.updatePhi;
            inputnames = {VarName({'..'}, 'phi')};
            fnmodel = {'..'};
            model = model.addPropFunction({'am', 'phiElectrode'}, fn, inputnames, fnmodel);
            
            fn = @ElectrodeActiveComponent.updateTemperature;
            fnmodel = {'..'};
            inputnames = {VarName(fnmodel, 'T')};
            model = model.addPropFunction({'am', 'T'}, fn, inputnames, fnmodel);
            
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
