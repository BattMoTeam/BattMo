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
            
            fn = @ProtonicMembraneCathode.updateJHp;
            inputnames = {'eta'};
            model = model.registerPropFunction({'jHp', fn, inputnames});


            % fn = @ProtonicMembraneCathode.updateJ;
            % inputnames = {'eta'};
            % model = model.registerPropFunction({'j', fn, inputnames});

            
            % fn = @ProtonicMembraneCathode.updatePi;
            % inputnames = {};
            % model = model.registerPropFunction({'pi', fn, inputnames});

            % model = model.removeVarName('eta');
            
        end
        
        
        function state = updateJHp(model, state)

            R = model.R;

            state.jHp = 1/R * state.eta;
            
        end


        % function state = updateJ(model, state)

        %     state.j = 1/model.R * state.eta;
            
        % end


        
        % function state = updatePi(model, state)

        %     state.pi = 0;
            
        % end
        
    end
    
end
