classdef Electrode_ < CompositeModel

    methods

        function model = Electrode_(name)
            
            model = model@CompositeModel(name);
            
            names = {model.names{:}, ...
                     'T' };
            
            model.SubModels{1} = ElectrodeActiveComponent_('eac');
            model.SubModels{2} = CurrentCollector_('cc');

            % Add coupling between the two components
            inputnames = {VarName({'..', 'eac'}, 'phi'), ...
                          VarName({'..', 'cc'}, 'phi')};
            fn = @Electrode.updateCoupling;
            fnmodel = {'..'};
            model = model.addPropFunction({'eac', 'jCoupling'}, fn, inputnames, fnmodel);
            model = model.addPropFunction({'cc', 'jCoupling'}, fn, inputnames, fnmodel);
            model = model.addPropFunction({'cc', 'eSource'}, fn, inputnames, fnmodel);
            
            % Temperature coupling between current collector and electrode active component
            inputnames = {VarName({'..', 'eac'}, 'T'), VarName({'..', 'cc'}, 'T')};
            fn = @Electrode.updateTemperatureCoupling;
            fnmodel = {'..'};
            model = model.addPropFunction({'eac', 'jHeatBcSource'}, fn, inputnames, fnmodel);
            model = model.addPropFunction({'cc', 'jHeatBcSource'}, fn, inputnames, fnmodel);
            
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
