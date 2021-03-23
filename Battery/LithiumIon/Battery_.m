classdef Battery_ < CompositeModel

    methods

        function model = Battery_()

            model = model@CompositeModel('battery');
            
            names = {'T', 'SOC'};
            model.names = names;
            
            submodels = {};
            submodels{end + 1} = orgLiPF6_('elyte');
            submodels{end + 1} = GraphiteElectrode_('ne');
            submodels{end + 1} = NMC111Electrode_('pe');
            submodels{end + 1} = CurrentCollector_('ccne');
            submodels{end + 1} = CurrentCollector_('ccpe');
            
            model.SubModels = submodels;
            model.hasparent = false;

            ccpe = model.getAssocModel('ccpe');
            ccpe.pnames = {ccpe.pnames{:}, 'E'};
            ccpe.names = {ccpe.names{:}, 'E'};
            ccpe.vardims('E') = 1;
            model = model.setSubModel('ccpe', ccpe);


            %% Setup property update functions 

            ne = model.getAssocModel('ne');
            pe = model.getAssocModel('pe');
            ccne = model.getAssocModel('ccne');
            ccpe = model.getAssocModel('ccpe');
            elyte = model.getAssocModel('elyte');

            % update function for temperature and soc
            fnupdate = @(model, state) model.dispatchValues(state);
            inputnames = {VarName({'..'}, 'T'), ...
                          VarName({'..'}, 'SOC')};
            fnmodel = {'..'};
            ne = ne.addPropFunction('T'  , fnupdate, inputnames, fnmodel);
            ne = ne.addPropFunction('SOC', fnupdate, inputnames, fnmodel);
            pe = pe.addPropFunction('T'  , fnupdate, inputnames, fnmodel);
            pe = pe.addPropFunction('SOC', fnupdate, inputnames, fnmodel);
            ccne = ccne.addPropFunction('T', fnupdate, inputnames, fnmodel);
            ccpe = ccpe.addPropFunction('T', fnupdate, inputnames, fnmodel);
            elyte = elyte.addPropFunction('T', fnupdate, inputnames, fnmodel);
            
            % update function for exchange terms (ne-elyte)
            fnupdate = @(model, state) setupExchanges(model, state);
            inputnames = {VarName({'..', 'ne'}, 'R'), ...
                          VarName({'..', 'pe'}, 'R')};
            fnmodel = {'..'};
            ne = ne.addPropFunction('LiSource', fnupdate, inputnames, fnmodel);
            ne = ne.addPropFunction('eSource' , fnupdate, inputnames, fnmodel);
            pe = pe.addPropFunction('LiSource', fnupdate, inputnames, fnmodel);
            pe = pe.addPropFunction('eSource' , fnupdate, inputnames, fnmodel);
            elyte = elyte.addPropFunction('LiSource', fnupdate, inputnames, fnmodel);
            
            % update function for phielyte (electrolyte potential)
            fnupdate = @(model, state) model.updatePhiElyte(state);
            inputnames = {VarName({'..', 'elyte'}, 'phi')};
            fnmodel = {'..'};
            ne = ne.addPropFunction('phielyte', fnupdate, inputnames, fnmodel);
            pe = pe.addPropFunction('phielyte', fnupdate, inputnames, fnmodel);
            
            % update function for boundary terms (ne-ccne)
            fnupdate = @(model, state) setupBCSources(model, state);
            inputnames = {VarName({'..', 'ne', 'am'}, 'phi'), ...
                          VarName({'..', 'pe', 'am'}, 'phi'), ... 
                          VarName({'..', 'ccne'}, 'phi'), ...
                          VarName({'..', 'ccpe'}, 'phi'), ...
                          VarName({'..', 'ccpe'}, 'E'), ...
                         };
            fnmodel = {'..'};
            ne   = ne.addPropFunction('jBcSource', fnupdate, inputnames, fnmodel);
            pe   = pe.addPropFunction('jBcSource', fnupdate, inputnames, fnmodel);
            ccne = ccne.addPropFunction('jBcSource', fnupdate, inputnames, fnmodel);
            ccpe = ccpe.addPropFunction('jBcSource', fnupdate, inputnames, fnmodel);
            
            
            fnupdate = @(model, state) model.dynamicBuildSOE(state); 
            % function above is not correct. This is just used now to inform the graph
            inputnames = {VarName({'..', 'ne'}   , 'LiSource'), ...
                          VarName({'..', 'pe'}   , 'LiSource'), ...
                          VarName({'..', 'elyte'}, 'LiSource'), ...
                          VarName({'..', 'ne'}   , 'LiFlux'), ...
                          VarName({'..', 'pe'}   , 'LiFlux'), ...
                          VarName({'..', 'elyte'}, 'LiFlux'), ...
                         };
            fnmodel = {'..'};
            ne    = ne.addPropFunction('massCont', fnupdate, inputnames, fnmodel);
            pe    = pe.addPropFunction('massCont', fnupdate, inputnames, fnmodel);
            elyte = elyte.addPropFunction('massCont', fnupdate, inputnames, fnmodel);
            
            model = model.setSubModel('ne', ne);
            model = model.setSubModel('pe', pe);
            model = model.setSubModel('ccne', ccne);
            model = model.setSubModel('ccne', ccne);
            model = model.setSubModel('elyte', elyte);

            model = model.initiateCompositeModel();
        end
        
    end
    
end
