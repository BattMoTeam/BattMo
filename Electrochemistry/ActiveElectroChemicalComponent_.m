classdef ActiveElectroChemicalComponent_ < ElectroChemicalComponent_

    methods

        function model = ActiveElectroChemicalComponent_(name)

            model = model@ElectroChemicalComponent_(name);
            
            model.SubModels{1} = ActiveMaterial_('am');
            
            propfunctions = model.propfunctions;
            
            fn = @Eletrode.updateIonSource;
            inputnames = {VarName({'am'}, 'R')};
            fnmodel = {'.'};
            propfunctions{end + 1} = PropFunction('chargeCarrierSource', fn, inputnames, fnmodel);
                        
            model.propfunctions = propfunctions;
            
            am = model.getAssocModel('am');
            
            fn = @ActiveElectroChemicalComponent.updateChargeCarrier;
            inputnames = {VarName({'..'}, 'cs')};
            fnmodel = {'..'};
            am = am.addPropFunction('cs', fn, inputnames, fnmodel);

            fn = @ActiveElectroChemicalComponent.updatePhi;
            inputnames = {VarName({'..'}, 'phi')};
            fnmodel = {'..'};
            am = am.addPropFunction('phi', fn, inputnames, fnmodel);

            fn = @ActiveElectroChemicalComponent.updateT;
            inputnames = {VarName({'..'}, 'T')};
            fnmodel = {'..'};
            am = am.addPropFunction('T', fn, inputnames, fnmodel);
            
            model = model.setSubModel('am', am);
            
        end
        
    end
    
end

