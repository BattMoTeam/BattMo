classdef ElectrodeThermalActiveComponent_ < ElectrodeActiveComponent_

    methods

        function model = ElectrodeActiveComponent_(name)

            model = model@ThermoElectroChemicalComponent_(name);
            
            names = {model.names{:}, ...
                     'jCoupling'};
            model.names = names;

            model.SubModels{1} = ActiveMaterial_('am');
            
            fn = @ElectrodeActiveComponent.updatejBcSource;
            inputnames = {'jCoupling'};
            fnmodel = {'.'};
            model = model.addPropFunction('jBcSource', fn, inputnames, fnmodel);            
            
            fn = @ElectrodeActiveComponent.updateIonAndCurrentSource;
            inputnames = {VarName({'am'}, 'R')};
            fnmodel = {'.'};
            model = model.addPropFunction('chargeCarrierSource', fn, inputnames, fnmodel);
            model = model.addPropFunction('eSource', fn, inputnames, fnmodel);
            
            fn = @ElectrodeActiveComponent.updateChargeCarrier;
            inputnames = {VarName({'..'}, 'chargeCarrier')};
            fnmodel = {'..'};
            model = model.addPropFunction({'am', 'cElectrode'}, fn, inputnames, fnmodel);

            fn = @ElectrodeActiveComponent.updatePhi;
            inputnames = {VarName({'..'}, 'phi')};
            fnmodel = {'..'};
            model = model.addPropFunction({'am', 'phiElectrode'}, fn, inputnames, fnmodel);
            
            fn = @ElectrodeActiveComponent.updateT;
            inputnames = {VarName({'..'}, 'T')};
            fnmodel = {'..'};
            model = model.addPropFunction({'am', 'T'}, fn, inputnames, fnmodel);
            
            fn = @ElectrodeActiveComponent.updateDiffusionCoefficient;
            inputnames = {VarName({'am'}, 'D')};
            fnmodel = {'.'};
            model = model.addPropFunction('D', fn, inputnames, fnmodel);
            
        end
        
    end
    
end

