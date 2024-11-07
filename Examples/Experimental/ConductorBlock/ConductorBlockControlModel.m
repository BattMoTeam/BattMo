classdef ConductorBlockControlModel < ControlModel

    properties

        Imax
        
    end
    
    methods

        function model = ConductorBlockControlModel(inputparams)

            model = model@ControlModel(inputparams);
            
            fdnames = {'Imax'};
            model = dispatchParams(model, inputparams, fdnames);
            
        end

        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@ControlModel(model);
            
            fn = @CCDischargeControlModel.updateControlEquation;
            model = model.registerPropFunction({'controlEquation', fn, {'I'}});
            
        end

        function state = updateControlEquation(model, state)
            
            Imax = model.Imax;
            I    = state.I;            

            state.controlEquation = I - Imax;
            
        end
        

    end
    
end
