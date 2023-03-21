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
        
        function paramobj = EvolutionElectrodeInputParams(jsonstruct)

            paramobj = paramobj@ComponentInputParams(jsonstruct);
            
            ptl = 'PorousTransportLayer';
            ctl = 'CatalystLayer';
            exr = 'ExchangeReaction';
            
            pick = @(fd) pickField(jsonstruct, fd);

            switch paramobj.porousTransportLayerType
              case {'Hydrogen', 'Oxygen'}
                % for the moment no difference in the input for both models
                paramobj.(ptl) = PorousTransportLayerInputParams(pick(ptl));
              otherwise
                error('porousTransportLayerType not recognized')
            end

            switch paramobj.catalystLayerType
              case {'Platinium', 'Iridium'}
                % for the moment no difference in the input for both models
                paramobj.(ctl) = CatalystLayerInputParams(pick(ctl));
              otherwise
                error('catalystLayerType not recognized')
            end
            
            paramobj.(exr) = ExchangeReactionInputParams(pick(exr));
            
        end
        
    end
    
    
end
