classdef Battery_ < CompositeModel

    methods

        function model = Battery_()

            model = model@CompositeModel('battery');
            
            names = {'T', ...
                     'SOC', ...
                     'prevstate', ...
                     'dt'};
            
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


            %% Declare property update functions 

            ne = model.getAssocModel('ne');
            pe = model.getAssocModel('pe');
            ccne = model.getAssocModel('ccne');
            ccpe = model.getAssocModel('ccpe');
            elyte = model.getAssocModel('elyte');

            %% Declare update function for temperature and soc
            fn = @Battery.dispatchValues;
            inputnames = {VarName({'..'}, 'T'), ...
                          VarName({'..'}, 'SOC')};
            fnmodel = {'..'};
            ne = ne.addPropFunction('T'  , fn, inputnames, fnmodel);
            ne = ne.addPropFunction('SOC', fn, inputnames, fnmodel);
            pe = pe.addPropFunction('T'  , fn, inputnames, fnmodel);
            pe = pe.addPropFunction('SOC', fn, inputnames, fnmodel);
            ccne = ccne.addPropFunction('T', fn, inputnames, fnmodel);
            ccpe = ccpe.addPropFunction('T', fn, inputnames, fnmodel);
            elyte = elyte.addPropFunction('T', fn, inputnames, fnmodel);

            inputnames = {VarName({'..', '..'}, 'T'), ...
                          VarName({'..', '..'}, 'SOC')};
            fnmodel = {'..', '..'};

            am = ne.getAssocModel('am');
            am = am.addPropFunction('T', fn, inputnames, fnmodel);
            am = am.addPropFunction('SOC', fn, inputnames, fnmodel);
            ne = ne.setSubModel('am', am);

            am = pe.getAssocModel('am');
            am = am.addPropFunction('T', fn, inputnames, fnmodel);
            am = am.addPropFunction('SOC', fn, inputnames, fnmodel);
            pe = pe.setSubModel('am', am);
            
            %% Declare update functions for accumulation term
            fn = @Battery.updatAccumTerm;
            inputnames = {'Li', {'..', 'prevstate'}, {'..', 'dt'}};
            fnmodel = {'.'};

            ne = ne.addPropFunction('LiAccum', fn, inputnames, fnmodel);
            pe = pe.addPropFunction('LiAccum', fn, inputnames, fnmodel);
            elyte = elyte.addPropFunction('LiAccum', fn, inputnames, fnmodel);
            
            %%  Declare update function for exchange terms (ne-elyte, pe-elyte)
            
            fn =  @Battery.setupExchanges;
            inputnames = {VarName({'..', 'ne'}, 'R'), ...
                          VarName({'..', 'pe'}, 'R')};
            fnmodel = {'..'};
            ne = ne.addPropFunction('LiSource', fn, inputnames, fnmodel);
            ne = ne.addPropFunction('eSource' , fn, inputnames, fnmodel);
            pe = pe.addPropFunction('LiSource', fn, inputnames, fnmodel);
            pe = pe.addPropFunction('eSource' , fn, inputnames, fnmodel);
            elyte = elyte.addPropFunction('LiSource', fn, inputnames, fnmodel);
            elyte = elyte.addPropFunction('eSource', fn, inputnames, fnmodel);
            ccne = ccne.addPropFunction('eSource', fn, inputnames, fnmodel);
            ccpe = ccpe.addPropFunction('eSource', fn, inputnames, fnmodel);
            
            %% Declare update function for phielyte (electrolyte potential is property of the electrode)
            
            fn =  @Battery.updatePhiElyte;
            inputnames = {VarName({'..', 'elyte'}, 'phi')};
            fnmodel = {'..'};
            ne = ne.addPropFunction('phielyte', fn, inputnames, fnmodel);
            pe = pe.addPropFunction('phielyte', fn, inputnames, fnmodel);
            
            %% Declare update functions for boundary terms (for ne, ccne, pe, ccpe)
            fn =  @Battery.setupBCSources;
            inputnames = {VarName({'..', 'ne', 'am'}, 'phi'), ...
                          VarName({'..', 'pe', 'am'}, 'phi'), ... 
                          VarName({'..', 'ccne'}, 'phi'), ...
                          VarName({'..', 'ccpe'}, 'phi'), ...
                          VarName({'..', 'ccpe'}, 'E'), ...
                         };
            fnmodel = {'..'};
            ne   = ne.addPropFunction('jBcSource', fn, inputnames, fnmodel);
            pe   = pe.addPropFunction('jBcSource', fn, inputnames, fnmodel);
            ccne = ccne.addPropFunction('jBcSource', fn, inputnames, fnmodel);
            ccpe = ccpe.addPropFunction('jBcSource', fn, inputnames, fnmodel);
            elyte = elyte.addPropFunction('jBcSource', fn, inputnames, fnmodel);
            
            model = model.setSubModel('ne', ne);
            model = model.setSubModel('pe', pe);
            model = model.setSubModel('ccne', ccne);
            model = model.setSubModel('ccpe', ccpe);
            model = model.setSubModel('elyte', elyte);

            model = model.initiateCompositeModel();
        end
        
    end
    
end
