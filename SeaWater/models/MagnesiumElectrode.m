classdef MagnesiumElectrode < SeaWaterElectrode

    methods
        
        function model = MagnesiumElectrode(inputparams)
            
            model = model@SeaWaterElectrode(inputparams);
            
        end
        
        function model = registerVarAndPropfuncNames(model)
            
            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)
            
            model = registerVarAndPropfuncNames@SeaWaterElectrode(model);
            
            model = model.registerVarName('E');
            
            fn = @() SeaWaterElectrode.updatePotential;
            inputnames = {'E'};
            model = model.registerPropFunction({'phi', fn, inputnames});
            
       end

       
        function state = updatePotential(model, state)
        % this is a special setup : We impose given potential and assume infinite conductivity in hydrogen electrode
        % (likely to change later)
            nc = model.G.cells.num;
            
            E = state.E;
            
            state.phi = E.*ones(nc, 1);
        end
        
        
    end
    
end