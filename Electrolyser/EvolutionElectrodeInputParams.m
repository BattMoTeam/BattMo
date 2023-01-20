classdef EvolutionElectrodeInputParams < InputParams
    
    properties
        
        PorousTransportLayer
        CatalystLayer
        ExchangeLayer

        porousTransportLayerType
        catalystLayerType

        couplingTerm
        
    end
    
    methods
        
        function paramobj = EvolutionElectrodeInputParams(jsonstruct)

            paramobj = paramobj@InputParams(jsonstruct);
            
            ptl = 'PorousTransportLayer';
            ctl = 'CatalystLayer';
            exl = 'ExchangeLayer';
            
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
            
            paramobj.(exl) = ExchangeLayerInputParams(pick(exl));
            
        end
        
    end
    
    
end
