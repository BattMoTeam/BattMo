classdef Scanning < handle

    properties (SetAccess = immutable)
        
        jsonstruct
        packingMass
        
    end

    properties

        model
        specificEnergyHour % in Wh/hour
        specificEnergy % 
        thickness 
        SiContent

    end
    
    methods

        function sc = Scanning()

            sc.packingMass = 25*gram;
            
            jsonstruct = parseBattmoJson('Examples/LNMO/SiGrLnmo.json');
            jsonstruct_geometry = parseBattmoJson('Examples/JsonDataFiles/geometryMultiLayerPouch.json');
            jsonstruct_geometry.Geometry.nLayers = 50;
            
            jsonstruct = mergeStructs({jsonstruct_geometry, ...
                                           jsonstruct});

            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            sep = 'Separator';
            cc  = 'CurrentCollector';
            co  = 'Coating';
            am1 = 'ActiveMaterial1';
            am2 = 'ActiveMaterial2';
            am  = 'ActiveMaterial';
            itf = 'Interface';

            % set width and length of pouch

            jsonstruct.Geometry.length = 300*milli*meter;
            jsonstruct.Geometry.width  = 90*milli*meter;
            
            % We change the guest stoichiometries
            
            jsonstruct.(ne).(co).(am1).(itf).guestStoichiometry100 = 0.99;
            jsonstruct.(ne).(co).(am1).(itf).guestStoichiometry0   = 0.01;
            jsonstruct.(ne).(co).(am2).(itf).guestStoichiometry100 = 0.99;
            jsonstruct.(ne).(co).(am2).(itf).guestStoichiometry0   = 0.01;

            jsonstruct.(pe).(co).(am).(itf).guestStoichiometry0   = 0.99;
            jsonstruct.(pe).(co).(am).(itf).guestStoichiometry100 = 0.01;

            % We do not need fine discretization, set it to minimum

            jsonstruct.(ne).(co).numberOfDiscreteCells = 1;
            jsonstruct.(pe).(co).numberOfDiscreteCells = 1;
            jsonstruct.(ne).(cc).numberOfDiscreteCells = 1;
            jsonstruct.(pe).(cc).numberOfDiscreteCells = 1;
            jsonstruct.(sep).numberOfDiscreteCells = 1;
            
            sc.jsonstruct = jsonstruct;
            
        end

        function computeSpecificEnergy(sc, thickness, SiContent)

            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            sep = 'Separator';
            co  = 'Coating';
            bd  = 'Binder';
            ad  = 'ConductingAdditive';
            am1 = 'ActiveMaterial1';
            am2 = 'ActiveMaterial2';
            am  = 'ActiveMaterial';
            itf = 'Interface';

            jsonstruct = sc.jsonstruct;
            
            jsonstruct.(ne).(co).thickness = thickness;
            jsonstruct.(ne).(co).(am2).massFraction = SiContent;
            jsonstruct.(ne).(co).(am1).massFraction = 1 - (jsonstruct.(ne).(co).(am2).massFraction + ...
                                                           jsonstruct.(ne).(co).(bd).massFraction + ...
                                                           jsonstruct.(ne).(co).(ad).massFraction);
            
            model = setupModelFromJson(jsonstruct);

            [capacity, capacities] = computeCellCapacity(model);

            npratio = capacities.(ne)/capacities.(pe);

            jsonstruct.(pe).(co).thickness = jsonstruct.(pe).(co).thickness*(npratio/1.1);

            model = setupModelFromJson(jsonstruct);
            
            energy = computeCellEnergy(model);
            mass   = computeCellMass(model, 'packingMass', sc.packingMass);

            specificEnergy = energy/mass;

            sc.specificEnergy     = specificEnergy;
            sc.specificEnergyHour = specificEnergy/hour;
            sc.thickness          = thickness;
            sc.SiContent          = SiContent;
            sc.model              = model;
        end
        
    end

    
end

%{
Copyright 2021-2026 SINTEF Industry, Sustainable Energy Technology
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
