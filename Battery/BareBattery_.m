classdef BareBattery_ < CompositeModel

    methods

        function model = BareBattery_()

            model = model@CompositeModel('battery');
            
            names = {'SOC', ...
                     'T', ...
                     'controlEq'};
            
            model.names = names;
            
            % setup the submodels
            submodels = {};
            % electrolyte
            submodels{end + 1} = Electrolyte_('elyte');
            % negative electrode
            ne_model = ElectrodeActiveComponent_('ne');
            ne_model.names{end + 1} = 'jExternal';
            submodels{end + 1} = ne_model;
            % positive electrode            
            pe_model = ElectrodeActiveComponent_('pe');
            pe_model.names{end + 1} = 'E';
            pe_model.names{end + 1} = 'I';
            pe_model.names{end + 1} = 'jExternal';
            submodels{end + 1} = pe_model;
            
            model.SubModels = submodels;
            
            %% update temperatures (dispatching)
            fn = @Battery.updateTemperature;
            
            fnmodel = {'..'};
            inputnames = {VarName({'..'}, 'T')};
            model = model.addPropFunction({'ne', 'T'}, fn, inputnames, fnmodel);
            model = model.addPropFunction({'pe', 'T'}, fn, inputnames, fnmodel);
            model = model.addPropFunction({'elyte', 'T'}, fn, inputnames, fnmodel);
           
                  
            %% setup couplings
            
            fn = @Battery.updateElectrodeCoupling;
            
            clear inputnames;
            inputnames{1} = VarName({'elyte'}, 'c');
            inputnames{1}.isNamingRelative = false;
            inputnames{2} = VarName({'elyte'}, 'phi');
            inputnames{2}.isNamingRelative = false;
            
            fnmodel = {'..', '..'};
            model = model.addPropFunction({'ne', 'am', 'phiElectrolyte'}, fn, inputnames, fnmodel);
            model = model.addPropFunction({'ne', 'am', 'cElectrolyte'}, fn, inputnames, fnmodel);
            model = model.addPropFunction({'pe', 'am', 'phiElectrolyte'}, fn, inputnames, fnmodel);
            model = model.addPropFunction({'pe', 'am', 'cElectrolyte'}, fn, inputnames, fnmodel);
            
            fn = @Battery.updateElectrolyteCoupling;
            
            clear inputnames;
            inputnames{1} = VarName({'ne', 'am'}, 'R');
            inputnames{1}.isNamingRelative = false;
            inputnames{2} = VarName({'pe', 'am'}, 'R');
            inputnames{2}.isNamingRelative = false;
            
            fnmodel = {'..'};
            model = model.addPropFunction({'elyte', 'massSource'}, fn, inputnames, fnmodel);
            model = model.addPropFunction({'elyte', 'eSource'}, fn, inputnames, fnmodel);
            
            
            fn = @Battery.setupExternalCouplingNegativeElectrode;
            inputnames = {'phi'};
            fnmodel = {'..'};
            model = model.addPropFunction({'ne', 'jBcSource'}, fn, inputnames, fnmodel);
            
            fn = @Battery.setupExternalCouplingPositiveElectrode;
            inputnames = {'phi', 'E'};
            fnmodel = {'.'};
            model = model.addPropFunction({'pe', 'jBcSource'}, fn, inputnames, fnmodel);
            
            %% setup control equation
            fn = @Batter.setupEIEquation;
            fnmodel = {'.'};
            inputnames = {{'pe', 'E'}, ...
                          {'pe', 'I'}, ...
                          {'pe', 'phi'}, ...
                         };
            model = model.addPropFunction({'controlEq'}, fn, inputnames, fnmodel);

                                                    
            %% setup external coupling at positive and negative electrodes
            
            fn = @Battery.setupExternalCouplingNegativeElectrode;
            inputnames = {'phi'};
            fnmodel = {'..'};
                      
            model = model.addPropFunction({'ne', 'jExternal'}, fn, inputnames, fnmodel);
            
            fn = @Battery.setupExternalCouplingPositiveElectrode;
            inputnames = {'phi', 'E'};
            fnmodel = {'.'};
            model = model.addPropFunction({'pe', 'jExternal'}, fn, inputnames, fnmodel);
            
            model = model.initiateCompositeModel();
            
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
