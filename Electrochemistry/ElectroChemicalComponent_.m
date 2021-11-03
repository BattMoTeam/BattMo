classdef ElectroChemicalComponent_ < ElectronicComponent_

    methods

        function model = ElectroChemicalComponent_(name)

            model = model@ElectronicComponent_(name);
            
            names = {model.names{:}        , ...
                     'c'
                     'massSource' , ...
                     'massFlux'   , ...
                     'massAccum'  , ...
                     'massCons'};
            model.names = names;
            
            model.vardims('cs') = 2;

            fn = @ElectroChemicalComponent.updateMassFlux;
            inputnames = {'c'};
            fnmodel = {'.'};
            model = model.addPropFunction('massFlux', fn, inputnames, fnmodel);        
            
            fn = @ElectroChemicalComponent.updateMassConservation;
            inputnames = {'massFlux', 'massSource', 'massAccum'};
            fnmodel = {'.'};
            model = model.addPropFunction('massCons', fn, inputnames, fnmodel);

        end
        
    end
    
end

