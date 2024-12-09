classdef ProtonicMembraneCellWithGasSupplyInputParams < ComponentInputParams
    
    properties

        T
        
        Cell
        GasSupply
        
        couplingTerm

    end
    
    methods
        
        function inputparams = ProtonicMembraneCellWithGasSupplyInputParams(jsonstruct)

            inputparams = inputparams@ComponentInputParams(jsonstruct);
            
            pick = @(fd) pickField(jsonstruct, fd);

            inputparams.Cell      = ProtonicMembraneCellInputParams(pick('Cell'));
            inputparams.GasSupply = ProtonicMembraneGasSupplyInputParams(pick('GasSupply'));
            
            inputparams = inputparams.validateInputParams();
            
            
        end
        
    end
    
end
