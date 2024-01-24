classdef EvolutionElectrodeInputParams < ComponentInputParams
    
    properties
        
        PorousTransportLayer
        CatalystLayer
        ExchangeReaction

        porousTransportLayerType
        catalystLayerType

        couplingTerm
        
    end
    
    methods
        
        function inputparams = EvolutionElectrodeInputParams(jsonstruct)

            inputparams = inputparams@ComponentInputParams(jsonstruct);
            
            ptl = 'PorousTransportLayer';
            ctl = 'CatalystLayer';
            exr = 'ExchangeReaction';
            
            pick = @(fd) pickField(jsonstruct, fd);

            switch inputparams.porousTransportLayerType
              case {'Hydrogen', 'Oxygen'}
                % for the moment no difference in the input for both models
                inputparams.(ptl) = PorousTransportLayerInputParams(pick(ptl));
              otherwise
                error('porousTransportLayerType not recognized')
            end

            switch inputparams.catalystLayerType
              case {'Platinium', 'Iridium'}
                % for the moment no difference in the input for both models
                inputparams.(ctl) = CatalystLayerInputParams(pick(ctl));
              otherwise
                error('catalystLayerType not recognized')
            end
            
            inputparams.(exr) = ExchangeReactionInputParams(pick(exr));
            
        end
        
    end
    
    
end



%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
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
