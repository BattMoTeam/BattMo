classdef ProtonicMembraneStokesGasSupplyBc < ProtonicMembraneGasSupplyBc
    
    
    methods
        
        function model = ProtonicMembraneStokesGasSupplyBc(inputparams)
            
            model = model@ProtonicMembraneGasSupplyBc();
            
        end

        function model = registerVarAndPropfuncNames(model)
            
            model = registerVarAndPropfuncNames@ProtonicMembraneGasSupplyBc(model);

            nGas = model.nGas;
            


        end
        
    end

    
end
