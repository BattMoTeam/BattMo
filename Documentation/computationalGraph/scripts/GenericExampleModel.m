classdef GenericExampleModel < BaseModel

    methods
        
        function model = registerVarAndPropfuncNames(model)
            
            model = registerVarAndPropfuncNames@BaseModel(model);

            n = 3;
            p = 2;
            q = 3;
            
            varnames = {};
            varnames{end + 1} = VarName({}, 'x', n);
            varnames{end + 1} = VarName({}, 'y', p);
            varnames{end + 1} = VarName({}, 'f', q);
            
            model = model.registerVarNames(varnames);
            
            fn = @ThermalModel.updateY1;
            outputvarname = VarName({}, 'y', p, 1);
            inputvarnames = {VarName({}, 'x', n, 1 : 2), VarName({}, 'y', n, 2)};

            model = model.registerPropFunction({outputvarname, fn, inputvarnames});

            fn = @ThermalModel.updateY2;
            outputvarname = VarName({}, 'y', p, 2);
            inputvarnames = {VarName({}, 'x', n, 1 : 3)};
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});

            fn = @ThermalModel.updateF1;
            outputvarname = VarName({}, 'f', q, 1);
            inputvarnames = {VarName({}, 'y', p, 1 : 2)};
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});

            fn = @ThermalModel.updateF2;
            outputvarname = VarName({}, 'f', q, 2);
            inputvarnames = {VarName({}, 'y', p, 1 : 2), VarName({}, 'x', n, 2)};
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});

            fn = @ThermalModel.updateF3;
            outputvarname = VarName({}, 'f', q, 3);
            inputvarnames = {VarName({}, 'x', p, [1, 3])};
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});

            
        end

    end
end
