classdef Battery_ < CompositeModel

    methods

        function model = Battery_()

            model = model@CompositeModel('battery');
            
            names = {'SOC', ...
                     'energyCons', ...
                     'controlEq'};
            
            model.names = names;
            
            submodels = {};
            submodels{end + 1} = Electrolyte_('elyte');
            submodels{end + 1} = Electrode_('ne');
            pe_model = Electrode_('pe');
            pe_model.names{end + 1} = 'E';
            pe_model.names{end + 1} = 'I';
            submodels{end + 1} = pe_model;
            submodels{end + 1} = ThermalComponent_('thermal');
            
            model.SubModels = submodels;
            
            %% update temperatures (dispatching)
            fn = @Battery.updateTemperature;
            
            fnmodel = {'..', '..'};
            inputnames = {VarName({'..', '..', 'thermal'}, 'T')};
            model = model.addPropFunction({'ne', 'eac', 'T'}, fn, inputnames, fnmodel);
            model = model.addPropFunction({'ne', 'cc', 'T'}, fn, inputnames, fnmodel);
            model = model.addPropFunction({'pe', 'eac', 'T'}, fn, inputnames, fnmodel);
            model = model.addPropFunction({'pe', 'cc', 'T'}, fn, inputnames, fnmodel);            
            fnmodel = {'..'};
            inputnames = {VarName({'..', 'thermal'}, 'T')};            
            model = model.addPropFunction({'elyte', 'T'}, fn, inputnames, fnmodel);
                  
            %% setup couplings
            
            fn = @Battery.updateElectrodeCoupling;
            
            clear inputnames;
            inputnames{1} = VarName({'elyte'}, 'c');
            inputnames{1}.isNamingRelative = false;
            inputnames{2} = VarName({'elyte'}, 'phi');
            inputnames{2}.isNamingRelative = false;
            
            fnmodel = {'..', '..', '..'};
            model = model.addPropFunction({'ne', 'eac', 'am', 'phiElectrolyte'}, fn, inputnames, fnmodel);
            model = model.addPropFunction({'ne', 'eac', 'am', 'cElectrolyte'}, fn, inputnames, fnmodel);
            model = model.addPropFunction({'pe', 'eac', 'am', 'phiElectrolyte'}, fn, inputnames, fnmodel);
            model = model.addPropFunction({'pe', 'eac', 'am', 'cElectrolyte'}, fn, inputnames, fnmodel);
            
            fn = @Battery.updateElectrolyteCoupling;
            
            clear inputnames;
            inputnames{1} = VarName({'ne', 'eac', 'am'}, 'R');
            inputnames{1}.isNamingRelative = false;
            inputnames{2} = VarName({'pe', 'eac', 'am'}, 'R');
            inputnames{2}.isNamingRelative = false;
            
            fnmodel = {'..'};
            model = model.addPropFunction({'elyte', 'cSource'}, fn, inputnames, fnmodel);
            model = model.addPropFunction({'elyte', 'eSource'}, fn, inputnames, fnmodel);
            
            %% setup control equation
            fn = @Batter.setupEIEquation;
            fnmodel = {'.'};
            inputnames = {{'pe', 'cc', 'E'}, ...
                          {'pe', 'cc', 'I'}, ...
                          {'pe', 'cc', 'phi'}, ...
                         };
            model = model.addPropFunction({'controlEq'}, fn, inputnames, fnmodel);
            
            %% update Thermal accumulation terms
            
            fn = @Battery.updateThermalAccumTerms;
            fnmodel = {'..'};
            inputnames = {'T'};
            model = model.addPropFunction({'thermal', 'accumHeat'}, fn, inputnames, fnmodel);
            
            %% update Thermal Ohmic Terms
            
            fn = @Battery.updateThermalOhmicSourceTerms;
            fnmodel = {'..'};
            inputnames = {VarName({'..', 'elyte'}, 'j')     , ...
                          VarName({'..', 'ne', 'cc'}, 'j')  , ...
                          VarName({'..', 'ne', 'eac'}, 'j') , ...
                          VarName({'..', 'pe', 'cc'}, 'j')  , ...
                          VarName({'..', 'pe', 'eac'}, 'j')};
            model = model.addPropFunction({'thermal', 'jHeatOhmSource'}, fn, inputnames, fnmodel);
            model = model.addPropFunction({'thermal', 'jHeatBcSource'}, fn, inputnames, fnmodel);
            
            %% update Thermal Chemical Terms
            
            fn = @Battery.updateThermalChemicalSourceTerms;
            fnmodel = {'..'};
            inputnames = {VarName({'..', 'elyte'}, 'diffFlux'), ...
                          VarName({'..', 'elyte'}, 'D'), ...
                          VarName({'..', 'elyte'}, 'dmudcs')};
            model = model.addPropFunction({'thermal', 'jHeatChemicalSource'}, fn, inputnames, fnmodel);
                          
            %% update Thermal Chemical Terms
            
            fn = @Battery.updateThermalReactionSourceTerms;
            fnmodel = {'..'};
            inputnames = {VarName({'..', 'ne', 'eac', 'am'}, 'R'), ...
                          VarName({'..', 'ne', 'eac', 'am'}, 'eta'), ...
                          VarName({'..', 'pe', 'eac', 'am'}, 'R'), ...
                          VarName({'..', 'pe', 'eac', 'am'}, 'eta')};
            model = model.addPropFunction({'thermal', 'jHeatReactionSource'}, fn, inputnames, fnmodel);
                                                    
            %% setup external coupling at positive and negative electrodes
            
            fn = @Battery.setupExternalCouplingNegativeElectrode;
            inputnames = {'phi'};
            fnmodel = {'..', '..'};
                      
            model = model.addPropFunction({'ne', 'cc', 'jExternal'}, fn, inputnames, fnmodel);
            
            fn = @Battery.setupExternalCouplingPositiveElectrode;
            inputnames = {'phi', ...
                          'E'};
            fnmodel = {'.'};
            model = model.addPropFunction({'pe', 'cc', 'jExternal'}, fn, inputnames, fnmodel);
            
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
