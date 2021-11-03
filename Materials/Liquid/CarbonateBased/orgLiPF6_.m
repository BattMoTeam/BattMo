classdef orgLiPF6_ < Electrolyte_

    methods

        function model = orgLiPF6_(name)

            model = model@Electrolyte_(name);

            fn = @orgLiPF6.updateCurrent;
            inputnames = {'Li', 'T', 'phi'};
            fnmodel = {'.'};
            model = model.addPropFunction('j', fn, inputnames, fnmodel);
            
            fn = @orgLiPF6.updateChemicalCurrent;
            inputnames = {'cs', 'j', 'T'};
            fnmodel = {'.'};
            model = model.addPropFunction('massFlux', fn, inputnames, fnmodel);        
            
        end
        
    end

end
