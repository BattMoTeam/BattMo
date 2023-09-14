classdef OxideMembraneElectrode < BaseModel
    
    properties

        constants
        
        T     % Temperature
        Rct   % charge transfer resistance
        Eocp  % Open circuit voltage
        muEl0 % standard value of chemical electron potential

        pO2   % O2 pressure used to compute Eocp
        
    end
    
    methods
        
        function model = OxideMembraneElectrode(paramobj)

            model = model@BaseModel();

            fdnames = {'T'     , ...
                       'pO2'   , ...
                       'muElO2', ...
                       'Rct'};
            model = dispatchParams(model, paramobj, fdnames);

            model.constants = PhysicalConstants();

            c = model.constants;
            
            Eocp = c.R*model.T/c.F*log(model.pO2);

            model.Eocp = Eocp;
            
        end
        
        function model = registerVarAndPropfuncNames(model)
        
            model = registerVarAndPropfuncNames@BaseModel(model);
            
            varnames = {};
            
            % Electrode electrostatic potential
            varnames{end + 1} = 'phi';
            % Electrode electromotive potential
            varnames{end + 1} = 'pi';
            % Electrode electronic chemical potential
            varnames{end + 1} = 'E';
            % Electron concentration
            varnames{end + 1} = 'ce';
            % Equation for pi - phi
            varnames{end + 1} = 'currentVoltageEq';
            % Electronic flux at electrode (positif is flux leaving electrode)
            varnames{end + 1} = 'jEl';
            % ionix flux at electrode (positif is flux leaving electrode)
            varnames{end + 1} = 'jO2';
            % Charge flux at electrode (positif is flux leaving electrode)
            varnames{end + 1} = 'j';
            % Charge conservation equation (that is j - jEl - jHp = 0)
            varnames{end + 1} = 'chargeCons';
            % Definition equation for jEl
            varnames{end + 1} = 'jElEquation';
            % Definition equation for jO2
            varnames{end + 1} = 'jO2Equation';
            % Electronic concentration definition equation
            varnames{end + 1} = 'elConcEq';
            model = model.registerVarNames(varnames);
        

            fn = @OxideMembraneElectrode.updateChargeCons;
            inputnames = {'j', 'jEl', 'jHp'};
            model = model.registerPropFunction({'chargeCons', fn, inputnames});

            fn = @OxideMembraneElectrode.updateE;
            inputnames = {'phi', 'pi'};
            model = model.registerPropFunction({'E', fn, inputnames});
            
            fn = @OxideMembraneElectrode.updateCurrentVoltageEquation;
            inputnames = {'j', 'E'};
            model = model.registerPropFunction({'chargeCons', fn, inputnames});

            fn = @OxideMembraneElectrode.updateElConcEquation;
            inputnames = {'pi', 'ce', 'phi'};
            model = model.registerPropFunction({'elConcEq', fn, inputnames});
                    
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
