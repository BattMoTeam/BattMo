classdef ProtonicMembraneControl < BaseModel
    
    properties

        controlType
        controlValue
        
    end
    
    
    methods
        
        function model = ProtonicMembraneControl(paramobj)

            model = model@BaseModel();

            fdnames = {'controlType', ...
                       'controlValue'};
            model = dispatchParams(model, paramobj, fdnames);
            
        end

        function model = registerVarAndPropfuncNames(model)

            varnames = {};

            % Potential
            varnames{end + 1} = 'U';
            % Current
            varnames{end + 1} = 'I';            
            % Control Equation
            varnames{end + 1} = 'controlEquation';

            model = model.registerVarNames(varnames);
            
            fn = @ProtonicMembraneControl.setupControl;
            switch model.controlType
              case 'current'
                inputnames = {'I'};
              case 'voltage'
                inputnames = {'U'};
              otherwise
                error('controlType not recognized');
            end
            model = model.registerPropFunction({'controlEquation', fn, inputnames});
            
        end

        function state = setupControl(model, state)

            switch model.controlType
              case 'current'
                state.controlEquation = state.I - model.controlValue;
              case 'voltage'
                state.controlEquation = state.U - model.controlValue;
              otherwise
                error('control type not recognized');
            end
            
        end

        
        
    end
    
end
