classdef orgLiPF6_ < ElectrochemicalComponent_

    methods

        function model = orgLiPF6_(name)

            model = model@ElectrochemicalComponent_(name);

            names = {'cs'        , ...
                     'LiSource'  , ...
                     'LiFlux'    , ...
                     'LiAccum'   , ...
                     'massCons'};
            model.names = names;
            
            model.vardims('cs') = 2;

            model = model.setAlias({'Li', VarName({'.'}, 'cs', 2, 1)});
            
            fn = @orgLiPF6.updateCurrent;
            inputnames = {'Li', 'T', 'phi'};
            fnmodel = {'.'};
            model = model.addPropFunction('j', fn, inputnames, fnmodel);
            
            fn = @orgLiPF6.updateLithiumFlux;
            inputnames = {'Li', 'j', 'T'};
            fnmodel = {'.'};
            model = model.addPropFunction('LiFlux', fn, inputnames, fnmodel);        
            
            fn = @orgLiPF6.updateMassConservation;
            inputnames = {'LiFlux', 'LiSource', 'LiAccum'};
            fnmodel = {'.'};
            model = model.addPropFunction('massCons', fn, inputnames, fnmodel);

        end
        
    end

end
