classdef ProtonicMembraneCellInputParams < ComponentInputParams
    
    properties

        T
        
        Electrolyser
        GasSupply
        
        couplingTerm

    end
    
    methods
        
        function inputparams = ProtonicMembraneCellInputParams(jsonstruct)

            inputparams = inputparams@ComponentInputParams(jsonstruct);
            
            pick = @(fd) pickField(jsonstruct, fd);

            inputparams.Electrolyser      = ProtonicMembraneInputParams(pick('Electrolyser'));
            inputparams.GasSupply = ProtonicMembraneGasSupplyInputParams(pick('GasSupply'));
            
            inputparams = inputparams.validateInputParams();
            
            
        end
        
    end
    
end
