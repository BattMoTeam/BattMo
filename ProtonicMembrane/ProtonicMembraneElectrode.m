classdef ProtonicMembraneElectrode < BaseModel
    
    properties
        
        T    % Temperature
        Eocp % Open circuit potential
        
        constants

        
    end
    
    methods
        
        function model = ProtonicMembraneElectrode(paramobj)

            model = model@BaseModel();

            fdnames = {'T', ...
                       'Eocp'};
            model = dispatchParams(model, paramobj, fdnames);
            
            model.constants = PhysicalConstants();
            
        end
        
        function model = registerVarAndPropfuncNames(model)
        
            model = registerVarAndPropfuncNames@BaseModel(model);
            
            varnames = {};
            
            % H+ flux at electrode (positif is flux leaving electrode)
            varnames{end + 1} = 'jHp';
            % Electronic flux at electrode (positif is flux leaving electrode)
            varnames{end + 1} = 'jEl';
            % Charge flux at electrode (positif is flux leaving electrode)
            varnames{end + 1} = 'j';
            % Open circuit potential at electrode
            varnames{end + 1} = 'Eocp';
            % Electrode electrostatic potential
            varnames{end + 1} = 'phi';
            % Electrode electromotive potential
            varnames{end + 1} = 'pi';
            % over potential
            varnames{end + 1} = 'eta';
            % Charge conservation equation (that is j - jEl - jHp = 0)
            varnames{end + 1} = 'chargeCons';
            % Definition equation for jEl
            varnames{end + 1} = 'jElEquation';
            % Definition equation for jHp
            varnames{end + 1} = 'jHpEquation';
            model = model.registerVarNames(varnames);
        

            fn = @ProtonicMembraneElectrode.updateChargeCons;
            inputnames = {'j', 'jEl', 'jHp'};
            model = model.registerPropFunction({'chargeCons', fn, inputnames});

            fn = @ProtonicMembraneElectrode.updateEta;
            inputnames = {'Eocp', 'phi', 'pi'};
            model = model.registerPropFunction({'eta', fn, inputnames});

            fn = @ProtonicMembraneElectrode.updateEocp;
            inputnames = {};
            model = model.registerPropFunction({'Eocp', fn, inputnames});            
                    
        end
        
        
        function state = updateChargeCons(model, state)

            j   = state.j;
            jEl = state.jEl;
            jHp = state.jHp;

            state.chargeCons = j - jEl - jHp;
           
        end

        function state = updateEta(model, state)

            Eocp = state.Eocp;
            pi   = state.pi;
            phi  = state.phi;
            
            state.eta = pi - phi - Eocp;
            
        end

        function state = updateEocp(model, state)

            state.Eocp = model.Eocp;
            
        end
        

    end
    
end
