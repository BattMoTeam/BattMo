classdef ProtonicMembraneCathode < ProtonicMembraneElectrode
    
    properties

        R_ct_H
        Ptot
        
    end
    
    methods
        
        function model = ProtonicMembraneCathode(paramobj)

            model = model@ProtonicMembraneElectrode(paramobj);

            fdnames = {'R_ct_H', ...
                       'Ptot'};
                       
            model = dispatchParams(model, paramobj, fdnames);

            c = model.constants;
            
            pH2O_neg = 0.05*model.Ptot; 
            pH2      = (model.Ptot - pH2O_neg)./2;
            
            model.Eocv = model.E_0 - (c.R*model.T/(2*c.F)).*log(pH2); % cathode half - cell potential
            
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

            R = model.R_ct_H;

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
