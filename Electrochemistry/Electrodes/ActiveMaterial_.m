classdef ActiveMaterial_ < CompositeModel
    
    methods
        
        function model = ActiveMaterial_(name)

            model = model@CompositeModel(name);
            names = {};
            
            % Temperature
            names{end + 1} = 'T';
            % Status of Charge
            names{end + 1} = 'SOC';
            % potential in electrode
            names{end + 1} = 'phiElectrode';
            % charge carrier concentration in electrode - value at surface
            names{end + 1} = 'cElectrode';
            % charge carrier concentration in electrode - Averaged value
            names{end + 1} = 'cElectrodeAveraged';
            % potential in electrolyte
            names{end + 1} = 'phiElectrolyte';
            % charge carrier concentration in electrolyte
            names{end + 1} = 'cElectrolyte';
            % eta
            names{end + 1} = 'eta';
            % Reaction rate
            names{end + 1} = 'R';
            % Diffusion coefficient
            % Note : This is the inner particle diffusion coefficient
            names{end + 1} = 'D';
            % OCP
            names{end + 1} = 'OCP';
            % Reaction rate coefficient
            names{end + 1} = 'j0';
            % Solid diffusion equation
            names{end + 1} = 'solidDiffusionEq';
            
            model.names = names;
            
            fn = @ActiveMaterial.updateReactionRateCoefficient;
            inputnames = {'cElectrode', 'cElectrolyte', 'T'};
            fnmodel = {'.'};
            model = model.addPropFunction('j0', fn, inputnames, fnmodel);

            fn = @ActiveMaterial.updateDiffusionCoefficient;
            inputnames = {'T'};
            fnmodel = {'.'};
            model = model.addPropFunction('D', fn, inputnames, fnmodel);

            fn = @ActiveMaterial.updateOCP;
            inputnames = {'cElectrode', 'T'};
            fnmodel = {'.'};
            model = model.addPropFunction('OCP', fn, inputnames, fnmodel);

            fn = @ActiveMaterial.updateReactionRate;
            inputnames = {'T', 'phiElectrolyte', 'phiElectrode', 'cElectrode', 'cElectrolyte', 'OCP', 'k'};
            fnmodel = {'.'};
            model = model.addPropFunction('R', fn, inputnames, fnmodel);
            % model = model.addPropFunction('eta', fn, inputnames, fnmodel);
            
            fn = @ActiveMaterial.assembleSolidDiffusionEquation;
            inputnames = {'D', 'R', 'cElectrode', 'cElectrodeAveraged'};
            fnmodel = {'.'};
            model = model.addPropFunction('solidDiffusionEq', fn, inputnames, fnmodel);
        
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
