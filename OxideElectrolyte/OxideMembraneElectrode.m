classdef OxideMembraneElectrode < BaseModel
    
    properties

        constants
        
        T     % Temperature
        Rct   % charge transfer resistance
        Eocp  % Open circuit voltage
        muEl0 % standard value of chemical electron potential
        Keh   % Equilibrium constant for the hole-electron reaction 
        
        pO2   % O2 pressure used to compute Eocp
        
    end
    
    methods
        
        function model = OxideMembraneElectrode(paramobj)

            model = model@BaseModel();

            fdnames = {'T'     , ...
                       'pO2'   , ...
                       'muElO2', ...
                       'Keh'   , ...
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
            % Hole concentration
            varnames{end + 1} = 'ch';
            % Equation for pi - phi
            varnames{end + 1} = 'currentVoltageEquation';
            % Equilibrium equation for hole-electron reaction
            varnames{end + 1} = 'equilibriumEquation';
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
            varnames{end + 1} = 'electronConcentrationEquation';
            model = model.registerVarNames(varnames);

            fn = @OxideMembraneElectrode.updateChargeCons;
            inputnames = {'j', 'jEl', 'jO2'};
            model = model.registerPropFunction({'chargeCons', fn, inputnames});

            fn = @OxideMembraneElectrode.updateE;
            inputnames = {'phi', 'pi'};
            model = model.registerPropFunction({'E', fn, inputnames});
            
            fn = @OxideMembraneElectrode.updateCurrentVoltageEquation;
            inputnames = {'j', 'E'};
            model = model.registerPropFunction({'currentVoltageEquation', fn, inputnames});

            fn = @OxideMembraneElectrode.updateElConcEquation;
            inputnames = {'pi', 'ce', 'phi'};
            model = model.registerPropFunction({'electronConcentrationEquation', fn, inputnames});

            fn = @OxideMembraneElectrode.updateEquilibriumReaction;
            inputnames = {'ch', 'ce'};
            model = model.registerPropFunction({'equilibriumEquation', fn, inputnames});

        end
        
        
        function state = updateChargeCons(model, state)

            j   = state.j;
            jEl = state.jEl;
            jO2 = state.jO2;

            state.chargeCons = j - jEl - jO2;
           
        end

        function state = updateE(model, state)

            pi   = state.pi;
            phi  = state.phi;
            
            state.E = pi - phi;
            
        end

        function state = updateCurrentVoltageEquation(model, state)

            Eocp = model.Eocp;
            Rct  = model.Eocp;

            E = state.E;
            j = state.j;
            
            cveq = E - Eocp - j*Rct;

            state.currentVoltageEquation = cveq;
            
        end

        function state = updateElConcEquation(model, state)

            mu0 = model.muEl0;
            c   = model.constants;
            
            ce = state.ce;
            E  = state.E;

            eceq = E + (1/c.F)*(mu0 + (c.R*c.T)*log(ce));

            state.electronConcentrationEquation = eceq;
            
        end

        
        function state = updateEquilibriumReaction(model, state)

            Keh = model.Keh;
            
            ch = state.ch;
            ce = state.ce;
            
            state.equilibriumEquation = ch*ce - Keh;
            
        end
    end
    
end
