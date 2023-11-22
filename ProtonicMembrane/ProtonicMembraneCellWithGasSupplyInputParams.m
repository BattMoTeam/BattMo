classdef ProtonicMembraneCellWithGasSupplyInputParams < ComponentInputParams
    
    properties

        T
        
        Cell
        GasSupply
        
        couplingTerms

    end
    
    methods
        
        function paramobj = ProtonicMembraneCellWithGasSupplyInputParams(jsonstruct)

            paramobj = paramobj@ComponentInputParams(jsonstruct);
            
            pick = @(fd) pickField(jsonstruct, fd);

            paramobj.Cell      = ProtonicMembraneCellInputParams(pick('Cell'));
            paramobj.GasSupply = ProtonicMembraneGasSupplyInputParams(pick('GasSupply'));
            
            paramobj = paramobj.validateInputParams();
            
            
        end
        
    end
    
end
