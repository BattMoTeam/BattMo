classdef ElectroChemicalComponent_ < ElectronicComponent_

    methods

        function model = ElectroChemicalComponent_(name)

            model = model@ElectronicComponent_(name);
            
            names = {model.names{:}        , ...
                     'cs'                  , ... % Species concentration (cell)
                     'D'                   , ... % Effective diffusion coefficient   
                     'chargeCarrierSource' , ...
                     'chargeCarrierFlux'   , ...
                     'chargeCarrierAccum'  , ...
                     'massCons'};
            model.names = names;
            
            model.vardims('cs') = 2;

            model = model.setAlias({'chargeCarrier', VarName({'.'}, 'cs', 2, 1)});
            
            fn = @ElectroChemicalComponent.updateDiffusionCoefficient;
            inputnames = {'T'};
            fnmodel = {'.'};
            model = model.addPropFunction('D', fn, inputnames, fnmodel);        

            fn = @ElectroChemicalComponent.updateChargeCarrierFlux;
            inputnames = {'D', 'chargeCarrier'};
            fnmodel = {'.'};
            model = model.addPropFunction('chargeCarrierFlux', fn, inputnames, fnmodel);        
            
            fn = @ElectroChemicalComponent.updateMassConservation;
            inputnames = {'chargeCarrierFlux', 'chargeCarrierSource', 'chargeCarrierAccum'};
            fnmodel = {'.'};
            model = model.addPropFunction('massCons', fn, inputnames, fnmodel);

        end
        
    end
    
end

