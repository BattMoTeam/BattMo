classdef ActiveMaterial_ < CompositeModel
    
    methods
        
        function model = ActiveMaterial_(name)

            model = model@CompositeModel(name);
            names = {};
            
            % Temperature
            names{end + 1} = 'T';
            % Status of Charge
            names{end + 1} = 'SOC';
            % potential in electrode
            names{end + 1} = 'phiElectrode';
            % charge carrier concentration in electrode
            % NOTE : this is value at the surface of the particle (as used in Butler-Volmer)
            names{end + 1} = 'cElectrode';
            % potential in electrolyte
            names{end + 1} = 'phiElectrolyte';
            % charge carrier concentration in electrolyte
            names{end + 1} = 'cElectrolyte';
            % eta
            names{end + 1} = 'eta';
            % Reaction rate
            names{end + 1} = 'R';
            % Diffusion coefficient
            % Note : This is the inner particle diffusion coefficient
            names{end + 1} = 'D';
            % OCP
            names{end + 1} = 'OCP';
            % Reaction constant
            names{end + 1} = 'k';
            
            model.names = names;
            
            % Alias: charge Carrier for cs{1}
            model = model.setAlias({'chargeCarrierElectrode', VarName({'.'}, 'cElectrode')});
            model = model.setAlias({'chargeCarrierElectrolyte', VarName({'.'}, 'cElectrolyte')});
            
            fn = @ActiveMaterial.updateDiffusionConductivityCoefficients;
            inputnames = {'T'};
            fnmodel = {'.'};
            model = model.addPropFunction('D', fn, inputnames, fnmodel);
            model = model.addPropFunction('k', fn, inputnames, fnmodel);
            
            fn = @ActiveMaterial.updateOCP;
            inputnames = {'chargeCarrierElectrode', 'T'};
            fnmodel = {'.'};
            model = model.addPropFunction('OCP', fn, inputnames, fnmodel);

            fn = @ActiveMaterial.updateReactionRate;
            inputnames = {'T', 'phiElectrolyte', 'phiElectrode', 'chargeCarrierElectrode', 'chargeCarrierElectrolyte', 'OCP', 'k'};
            fnmodel = {'.'};
            model = model.addPropFunction('R', fn, inputnames, fnmodel);
            model = model.addPropFunction('eta', fn, inputnames, fnmodel);
            
        end        
        
    end
    
end
