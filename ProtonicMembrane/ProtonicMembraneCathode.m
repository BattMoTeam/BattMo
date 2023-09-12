classdef ProtonicMembraneCathode < ProtonicMembraneElectrode
    
    properties

        R
        
    end
    
    methods
        
        function model = ProtonicMembraneCathode(paramobj)

            model = model@ProtonicMembraneElectrode(paramobj);

            fdnames = {'R'};
            model = dispatchParams(model, paramobj, fdnames);
            
        end
        
        function model = registerVarAndPropfuncNames(model)
        
            model = registerVarAndPropfuncNames@ProtonicMembraneElectrode(model);
            
            fn = @ProtonicMembraneCathode.updateJ;
            inputnames = {'eta'};
            model = model.registerPropFunction({'j', fn, inputnames});

                    
        end
        
        
        function state = updateJ(model, state)

            state.j = 1/model.R * state.eta;

        end


        
    end
    
end