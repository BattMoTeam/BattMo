classdef ProtonicMembraneElectrode < BaseModel
    
    properties
        
        T    % Temperature
        E_0  % Standard potential
        Eocv % Open circuit potential (value depends on conditions at electrode)

        N % discretization constant (number of values)
        
        constants

        gasSupplyType % Either 'coupled' or 'notCoupled'
        nGas          % Number of gas (each of them will have a partial pressure). Only needed when gasSupplyType == 'coupled'
        gasInd        % Structure whose fieldname give index number of the corresponding gas component.
        
    end
    
    methods
        
        function model = ProtonicMembraneElectrode(paramobj)

            model = model@BaseModel();

            fdnames = {'T'  , ...
                       'E_0', ...
                       'N'  , ...
                       'gasSupplyType'};
            model = dispatchParams(model, paramobj, fdnames);

            if isempty(model.gasSupplyType)
                model.gasSupplyType = 'notCoupled';
            end
            
            model.constants = PhysicalConstants();
            
        end
        
        function model = registerVarAndPropfuncNames(model)
        
            model = registerVarAndPropfuncNames@BaseModel(model);
            
            varnames = {};
            
            % H+ flux at electrode (positif is flux leaving electrode). Unit A (one value per cell)
            varnames{end + 1} = 'iHp';
            % Electronic flux at electrode (positif is flux leaving electrode) Unit A (one value per cell)
            varnames{end + 1} = 'iEl';
            % Charge flux at electrode (positif is flux leaving electrode) Unit A (one value per cell)
            varnames{end + 1} = 'i';
            % Open circuit potential at electrode
            varnames{end + 1} = 'Eocv';
            % Electrode electrostatic potential
            varnames{end + 1} = 'phi';
            % Electrode electromotive potential
            varnames{end + 1} = 'pi';
            % over potential
            varnames{end + 1} = 'eta';
            % Charge conservation equation (that is j - iEl - iHp = 0)
            varnames{end + 1} = 'chargeCons';
            % Definition equation for iEl
            varnames{end + 1} = 'iElEquation';
            % Definition equation for iHp
            varnames{end + 1} = 'iHpEquation';

            switch model.gasSupplyType
              case 'notCoupled'
                % nothing to add
              case 'coupled'
                varnames{end + 1} = VarName({}, 'pressures', model.nGas);
              otherwise
                error('gasSupplyType');
            end
            
            model = model.registerVarNames(varnames);
        

            fn = @ProtonicMembraneElectrode.updateChargeCons;
            inputnames = {'i', 'iEl', 'iHp'};
            model = model.registerPropFunction({'chargeCons', fn, inputnames});

            fn = @ProtonicMembraneElectrode.updateEta;
            inputnames = {'Eocv', 'phi', 'pi'};
            model = model.registerPropFunction({'eta', fn, inputnames});

            fn = @ProtonicMembraneElectrode.updateEocv;
            inputnames = {};
            model = model.registerPropFunction({'Eocv', fn, inputnames});            
                    
        end
        
        
        function state = updateChargeCons(model, state)

            i   = state.i;
            iEl = state.iEl;
            iHp = state.iHp;

            state.chargeCons = i - iEl - iHp;
           
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
