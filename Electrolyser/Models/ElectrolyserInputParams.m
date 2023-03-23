classdef ElectrolyserInputParams < InputParams
    
    properties

        G % parent grid (the component grids are subgrid of that one)

        Eref % Reference potential

        IonomerMembrane
        HydrogenEvolutionElectrode
        OxygenEvolutionElectrode        
                
        couplingTerms
        
        controlI % given value for galvanistic control
    end
    
    methods
        
        function paramobj = ElectrolyserInputParams(jsonstruct)

            paramobj = paramobj@InputParams(jsonstruct);
            
            pick = @(fd) pickField(jsonstruct, fd);

            paramobj.IonomerMembrane = IonomerMembraneInputParams(pick('IonomerMembrane'));
            paramobj.HydrogenEvolutionElectrode = EvolutionElectrodeInputParams(pick('HydrogenEvolutionElectrode'));
            paramobj.OxygenEvolutionElectrode = EvolutionElectrodeInputParams(pick('OxygenEvolutionElectrode'));

            paramobj = paramobj.validateInputParams();
        end

        function paramobj = validateInputParams(paramobj)

            assert(strcmp(paramobj.HydrogenEvolutionElectrode.porousTransportLayerType, 'Hydrogen'), 'Expected porous transport layer is hydrogen');
            assert(strcmp(paramobj.OxygenEvolutionElectrode.porousTransportLayerType, 'Oxygen'), 'Expected porous transport layer is oxygen');

            paramobj = mergeParameters(paramobj, {{'Eref'}                                               , ...
                                                  {'HydrogenEvolutionElectrode', 'CatalystLayer', 'Eref'}, ...
                                                  {'OxygenEvolutionElectrode', 'CatalystLayer', 'Eref'}});
            
            paramobj = validateInputParams@InputParams(paramobj);
            
        end
        
    end
    
    
end



%{
Copyright 2021-2023 SINTEF Industry, Sustainable Energy Technology
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
