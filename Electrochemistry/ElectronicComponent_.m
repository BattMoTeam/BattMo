classdef ElectronicComponent_ < CompositeModel

    methods

        function model = ElectronicComponent_(name)

            model = model@CompositeModel(name);
            names = {'T'        , ...
                     'phi'      , ...
                     'jBcSource', ...
                     'eSource'  , ...
                     'j'        , ...
                     'chargeCons'};
            model.names = names;
            
            propfunctions = model.propfunctions;
            
            fn = @ElectrochemicalComponent.updateCurrent;
            inputnames = {'phi'};
            fnmodel = {'.'};
            propfunctions{end + 1} = PropFunction('j', fn, inputnames, fnmodel);
            
            fn = @ElectrochemicalComponent.updateChargeConservation;
            inputnames = {'j', 'jBcSource', 'eSource'};
            fnmodel = {'.'};
            propfunctions{end + 1} = PropFunction('chargeCons', fn, inputnames, fnmodel);
            
            model.propfunctions = propfunctions;
            
        end
        
    end
    
end

       