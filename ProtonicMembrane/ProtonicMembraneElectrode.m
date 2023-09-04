classdef ProtonicMembraneElectrode < BaseModel
    
    properties
        
        T  % Temperature
        constants
        
        % coefficient in Buttler-Volmer
        beta
        % Exchange current density
        iBV_0
        % Limiting current densities
        anodeLimitingCurrentDensity
        cathodeLimitingCurrentDensity
        % charge
        z
        
    end
    
    methods
        
        function model = ProtonicMembraneElectrode(paramobj)

            model = model@BaseModel();

            fdnames = {'G', ...
                       'beta', ...
                       'iBV_0', ...
                       'anodeLimitingCurrentDensity', ...
                       'cathodeLimitingCurrentDensity', ...
                       'z'};
            
            % model.operators = localSetupOperators(model.G);

            model.constants = PhysicalConstants();
            
        end
        
        function model = registerVarAndPropfuncNames(model)
        
            model = registerVarAndPropfuncNames@BaseModel(model);
            
            varnames = {};
            
            % H+ flux at electrode
            varnames{end + 1} = 'jHp';
            % Electronic flux at electrode
            varnames{end + 1} = 'jEl';
            % Charge flux at electrode
            varnames{end + 1} = 'j';
            % Open circuit potential at electrode
            varnames{end + 1} = 'Eocp';
            % Electrode electrostatic potential
            varnames{end + 1} = 'phi';
            % Electrode electromotive potential
            varnames{end + 1} = 'pi';
            % over potential
            varnames{end + 1} = 'eta';
            % Charge conservation equation (j = jEl + jHp)
            varnames{end + 1} = 'chargeCons';
            % Definition equation for jEl
            varnames{end + 1} = 'jElEquation';
            % Definition equation for jHpl
            varnames{end + 1} = 'jHpEquation';
            model = model.registerVarNames(varnames);
        

            fn = @ProtonicMembraneElectrode.updateChargeCons;
            inputnames = {'j', 'jEl', 'jHp'};
            model = model.registerPropFunction({'chargeCons', fn, inputnames});

            fn = @ProtonicMembraneElectrode.updateEta;
            inputnames = {'Eocp', 'phi', 'pi'};
            model = model.registerPropFunction({'eta', fn, inputnames});
            
            fn = @ProtonicMembraneElectrode.updateJHp;
            inputnames = {'eta'};
            model = model.registerPropFunction({'jHp', fn, inputnames});

            model = model.registerStaticVarName('Eocp');
                    
        end
        
        
        function state = updateButtlerVolmerRate(model, state)
            
            R  = model.constants.R;
            F  = model.constants.F;
            T  = model.T;
            z  = model.z;
            ia = model.anodeLimitingCurrentDensity;
            ic = model.cathodeLimitingCurrentDensity;
            i0 = model.iBV_0;
            
            Eocp = state.Eocp
            E    = state.E;
            
            f = F*z/(RT);
            eta = E - Eocp
            
            iBV = i0*(exp(-beta*f*eta) - exp((1 - beta)*f*eta))/(1 + (i0/ic)*exp(-beta*f*eta) - (i0/ia)*exp(-(1 - beta)*f*eta));
            
            state.iBV = iBV;
            
        end
        
        
        
        
        
    end
    
end
