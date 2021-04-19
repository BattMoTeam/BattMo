classdef Electrolyte_ < ThermoElectroChemicalComponent_
%
% The function that couples the electolyte and the electrode is not implemented here (but at the battery level)
%
%
    methods
        
        function model = Electrolyte_(name)
            
            model = model@ThermoElectroChemicalComponent_(name);
            
            names = {model.names{:}, ...
                     'kappa'       , ...
                     'jchems'};
            model.names = names;
            
            model.vardims('dmucs') = 2;
            model.vardims('jchems') = 2;
            
            fn = @Electrolyte.updateChemicalCurrent;
            inputnames = {'chargeCarrier', 'T', 'phi'};
            fnmodel = {'.'};
            model = model.addPropFunction('kappa', fn, inputnames, fnmodel);
            model = model.addPropFunction('jchems', fn, inputnames, fnmodel);

            fn = @Electrolyte.updateDiffusionCoefficient;
            inputnames = {'chargeCarrier', 'T'};
            fnmodel = {'.'};
            model = model.addPropFunction('D', fn, inputnames, fnmodel);

            fn = @Electrolyte.updateCurrent;
            inputnames = {'phi', 'jchems', 'kappa'};
            fnmodel = {'.'};
            model = model.addPropFunction('j', fn, inputnames, fnmodel);
            
            fn = @Electrolyte.updateChargeCarrierFlux;
            inputnames = {'chargeCarrier', 'j', 'D'};
            fnmodel = {'.'};
            model = model.addPropFunction('chargeCarrierFlux', fn, inputnames, fnmodel);

            fn = @Electrolyte.updateCurrentBcSource;
            inputnames = {''};
            fnmodel = {'.'};
            model = model.addPropFunction('jBcSource', fn, inputnames, fnmodel);
            
        end        
    end
    
end