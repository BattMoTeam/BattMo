classdef ReactiveComponent_ < ElectronicComponent_

    methods

        function model = ReactiveComponent_(name)

            model = model@ElectroChemicalCompoent_(name);
            
            names = {model.names{:}, ...
                     'SOC'           , ...
                     'phielyte', ...
                     'R'             , ...
                     'LiSource'      , ...
                     'LiFlux'        , ...
                     'LiAccum'       , ...
                     'massCons'};
            model.names = names;
            
            model.SubModels{1} = ActiveComponent_('am');
            
            propfunctions = model.propfunctions;
            
            fn = @Electrode.updatePhi;
            inputnames = {VarName({'am'}, 'phi')};
            fnmodel = {'am'};
            propfunctions{end + 1} = PropFunction('phi', fn, inputnames, fnmodel);

            fn = @Eletrode.updateIonFlux;
            inputnames = {VarName({'am'}, 'D'), VarName({'am'}, 'Li')};
            fnmodel = {'.'};
            propfunctions{end + 1} = PropFunction('LiFlux', fn, inputnames, fnmodel);
            
            fn = @Electrode.updateMassConservation;
            inputnames = {'LiFlux', 'LiSource', 'LiAccum'};
            fnmodel = {'.'};
            propfunctions{end + 1} = PropFunction('massCons', fn, inputnames, fnmodel);
            
            fn = @Electrode.updateReactionRate;
            inputnames = {'T', 'phielyte', VarName({'am'}, 'phi'), VarName({'am'}, 'OCP'), VarName({'am'}, 'k')};
            fnmodel = {'.'};
            propfunctions{end + 1} = PropFunction('R', fn, inputnames, fnmodel);

            model.propfunctions = propfunctions;
        end
        
    end
    
end

