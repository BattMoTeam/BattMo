classdef ElectrochemicalComponent_ < ComponentModel

    methods

        function model = ElectrochemicalComponent_(name)

            model = model@ComponentModel(name);
            names = {'T'        , ...
                     'phi'      , ...
                     'jBcSource', ...
                     'eSource'  , ...
                     'j'        , ...
                     'chargeCons'};

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

       