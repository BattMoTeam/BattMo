classdef ProtonicMembraneCathode < ProtonicMembraneElectrode
    
    properties

        R_ct_H
        Ptot
        
        pRef

    end
    
    methods
        
        function model = ProtonicMembraneCathode(inputparams)

            model = model@ProtonicMembraneElectrode(inputparams);

            fdnames = {'R_ct_H', ...
                       'Ptot'};
                       
            model = dispatchParams(model, inputparams, fdnames);

            model.pRef = 1*barsa;

            c = model.constants;
            
            pH2O_neg = 0.1*model.Ptot; 
            pH2      = (model.Ptot - pH2O_neg);

            pRef = model.pRef;
            
            model.Eocv = model.E_0 - (c.R*model.T/(2*c.F)).*log(pH2/pRef); % cathode half - cell potential
            
        end
        
        function model = registerVarAndPropfuncNames(model)
        
            model = registerVarAndPropfuncNames@ProtonicMembraneElectrode(model);
            
            fn = @ProtonicMembraneCathode.updateIHp;
            inputnames = {'eta'};
            model = model.registerPropFunction({'iHp', fn, inputnames});


            % fn = @ProtonicMembraneCathode.updateJ;
            % inputnames = {'eta'};
            % model = model.registerPropFunction({'i', fn, inputnames});

            
            % fn = @ProtonicMembraneCathode.updatePi;
            % inputnames = {};
            % model = model.registerPropFunction({'pi', fn, inputnames});

            % model = model.removeVarName('eta');
            
        end
        
        
        function state = updateIHp(model, state)

            R = model.R_ct_H;

            state.iHp = 1/R * state.eta;
            
        end


        % function state = updateJ(model, state)

        %     state.i = 1/model.R * state.eta;
            
        % end


        
        % function state = updatePi(model, state)

        %     state.pi = 0;
            
        % end
        
    end
    
end
