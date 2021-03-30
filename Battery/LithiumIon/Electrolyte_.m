classdef Electrolyte_ < ElectroChemicalComponent_
%
% The function that couples the electolyte and the electrode is not implemented here (but at the battery level)
%
%
    function model = ElectroChemicalComponent_(name)
        
        names = {model.names{:}, ...
                 'kappa',
                 'dmudcs',
                 'jchems'};
        model.names = names;
        
        model.vardims('dmucs') = 2;
        model.vardims('jchems') = 2;
        
        fn = @Electrolyte.updateChemicalCurrent;
        inputnames = {'chargeCarrier', 'T', 'phi'};
        fnmodel = {'.'};
        model = model.addPropFunction('kappa', fn, inputnames, fnmodel);
        model = model.addPropFunction('dmudcs', fn, inputnames, fnmodel);
        model = model.addPropFunction('jchems', fn, inputnames, fnmodel);
    end
    
    
end