classdef GraphiteElectrode_ < Electrode_

    methods

        function model = GraphiteElectrode_(name)

            model = model@Electrode_(name);

            submodels = {};
            submodels{end + 1} = Graphite_('am');
            model.SubModels = submodels;

            model = model.setAlias({'Li', VarName({'am'}, 'Li')});

            am = model.getAssocModel('am');

            fn = @UpdateT;
            fnmodel = {'.'};
            inputvarnames = {VarName({'..'}, 'T')};
            am = am.addPropFunction('T', fn, inputvarnames, fnmodel);
            
            fn = @UpdateSOC;
            fnmodel = {'.'};
            inputvarnames = {VarName({'..'}, 'SOC')};
            am = am.addPropFunction('SOC', fn, inputvarnames, fnmodel);
            
            model = model.setSubModel('am', am);
           
        end
        
        
    end
    
end

       