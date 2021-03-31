classdef ActiveElectroChemicalComponent_ < ElectroChemicalComponent_

    methods

        function model = ActiveElectroChemicalComponent_(name)

            model = model@ElectroChemicalComponent_(name);
            
            names = {model.names{:}, ...
                     'jCoupling'};
            model.names = names;

            model.SubModels{1} = ActiveMaterial_('am');
            
            fn = @ActiveElectroChemicalComponent.updatejBcSource;
            inputnames = {'jCoupling'};
            fnmodel = {'.'};
            model = model.addPropFunction('jBcSource', fn, inputnames, fnmodel);            
            
            fn = @ActiveElectroChemicalComponent.updateIonAndCurrentSource;
            inputnames = {VarName({'am'}, 'R')};
            fnmodel = {'.'};
            model = model.addPropFunction('chargeCarrierSource', fn, inputnames, fnmodel);
            model = model.addPropFunction('eSource', fn, inputnames, fnmodel);
            
            fn = @ActiveElectroChemicalComponent.updateChargeCarrier;
            inputnames = {VarName({'am'}, 'cs')};
            fnmodel = {'.'};
            model = model.addPropFunction('cs', fn, inputnames, fnmodel);

            fn = @ActiveElectroChemicalComponent.updatePhi;
            inputnames = {VarName({'am'}, 'phi')};
            fnmodel = {'.'};
            model = model.addPropFunction('phi', fn, inputnames, fnmodel);
            
            fn = @ActiveElectroChemicalComponent.updateT;
            inputnames = {VarName({'..'}, 'T')};
            fnmodel = {'..'};
            model = model.addPropFunction({'am', 'T'}, fn, inputnames, fnmodel);
            
            fn = @ActiveElectroChemicalComponent.updateDiffusionCoefficient;
            inputnames = {VarName({'am'}, 'D')};
            fnmodel = {'.'};
            model = model.addPropFunction('D', fn, inputnames, fnmodel);
            
        end
        
    end
    
end

