classdef ProtonicMembraneElectrode < BaseModel
    
    properties
        
        T    % Temperature
        E_0  % Standard potential
        Eocv % Open circuit potential (value depends on conditions at electrode)
        
        constants
        
    end
    
    methods
        
        function model = ProtonicMembraneElectrode(paramobj)

            model = model@BaseModel();

            fdnames = {'T', ...
                       'E_0'};
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
            varnames{end + 1} = 'Eocv';
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
            inputnames = {'Eocv', 'phi', 'pi'};
            model = model.registerPropFunction({'eta', fn, inputnames});

            fn = @ProtonicMembraneElectrode.updateEocv;
            inputnames = {};
            model = model.registerPropFunction({'Eocv', fn, inputnames});            
                    
        end
        
        
        function state = updateChargeCons(model, state)

            j   = state.j;
            jEl = state.jEl;
            jHp = state.jHp;

            state.chargeCons = j - jEl - jHp;
           
        end

        function state = updateEta(model, state)

            Eocv = state.Eocv;
            pi   = state.pi;
            phi  = state.phi;
            
            state.eta = pi - phi - Eocv;
            
        end

        function state = updateEocv(model, state)

            state.Eocv = model.Eocv;
            
        end
        

    end
    
end
