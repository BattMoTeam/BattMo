classdef EvolutionElectrodeInputParams < InputParams
    
    properties
        
        PorousTransportLayer
        CatalystLayer
        ExchangeLayer

        porousTransportLayerType
        
    end
    
    methods
        
        function paramobj = EvolutionElectrodeInputParams(jsonstruct)

            paramobj = paramobj@InputParams(jsonstruct);
            
            ptl = 'PorousTransportLayer';
            ctl = 'CatalystLayer';
            exl = 'ExchangeLayer';
            
            pick = @(fd) pickField(jsonstruct, fd);

            switch paramobj.porousTransportLayerType
              case 'Hydrogen'
                paramobj.(ptl) = HydrogenPorousTransportLayerInputParams(pick(ptl));
              case 'Oxygen'
                paramobj.(ptl) = OxygenPorousTransportLayerInputParams(pick(ptl));
              otherwise
                error('porousTransportLayerType not recognized')
            end
            
            paramobj.(ctl) = CatalystLayerInputParams(pick(ctl));
            paramobj.(exl) = ExchangeLayerInputParams(pick(exl));
        end
        
    end
    
    
end
