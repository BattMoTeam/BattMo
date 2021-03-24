classdef Electrode_ < ElectrochemicalComponent_

    methods

        function model = Electrode_(name)

            model = model@ElectrochemicalComponent_(name);
            
            names = {model.names{:}, ...
                     'SOC'           , ...
                     'phiElectrolyte', ...
                     'R'             , ...
                     'LiSource'      , ...
                     'LiFlux'        , ...
                     'LiAccum'       , ...
                     'massCons'};
            model.names = names;
            
            propfunctions = model.propfunctions;
            
            fn = @updatePhi;
            inputnames = {VarName({'am'}, 'phi')};
            fnmodel = {'am'};
            propfunctions{end + 1} = PropFunction('phi', fn, inputnames, fnmodel);

            fn = @Eletrode.updateIonFlux;
            inputnames = {VarName({'am'}, 'D'), VarName({'am'}, 'cLi')};
            fnmodel = {'.'};
            propfunctions{end + 1} = PropFunction('LiFlux', fn, inputnames, fnmodel);
            
            fn = @Electrode.updateMassConservation;
            inputnames = {'LiFlux', 'LiSource', 'LiAccum'};
            fnmodel = {'.'};
            propfunctions{end + 1} = PropFunction('massCons', fn, inputnames, fnmodel);
            
            fn = @Electrode.updateReactionRate;
            inputnames = {'T', 'phiElectrolyte', VarName({'am'}, 'phi'), VarName({'am'}, 'OCP'), VarName({'am'}, 'k')};
            fnmodel = {'.'};
            propfunctions{end + 1} = PropFunction('R', fn, inputnames, fnmodel);

            model.propfunctions = propfunctions;
        end
        
    end
    
end

       