classdef OxideMembraneControl < BaseModel
    
    properties

        controlType
        
    end
    
    
    methods
        
        function model = OxideMembraneControl(inputparams)

            model = model@BaseModel();

            fdnames = {'controlType'};
            model = dispatchParams(model, inputparams, fdnames);
            
        end

        function model = registerVarAndPropfuncNames(model)

            varnames = {};

            % Potential
            varnames{end + 1} = 'U';
            % Current
            varnames{end + 1} = 'I';            
            % Value of the variable that is controled
            varnames{end + 1} = 'ctrlVal';
            % Control Equation
            varnames{end + 1} = 'controlEquation';

            model = model.registerVarNames(varnames);
            
            fn = @OxideMembraneControl.setupControl;
            switch model.controlType
              case 'current'
                inputnames = {'I', 'ctrlVal'};
              case 'voltage'
                inputnames = {'U', 'ctrlVal'};
              otherwise
                error('controlType not recognized');
            end
            model = model.registerPropFunction({'controlEquation', fn, inputnames});

        end

        function state = setupControl(model, state)

            switch model.controlType
              case 'current'
                state.controlEquation = state.I - state.ctrlVal;
              case 'voltage'
                state.controlEquation = state.U - state.ctrlVal;
              otherwise
                error('control type not recognized');
            end
            
        end
        
    end
    
end
