classdef MagnesiumAqueousModel < AqueousModel

    methods
        
        function model = MagnesiumAqueousModel(elyteModel, totals, totalnames, pH)

            model = model@AqueousModel(elyteModel, totals, totalnames, pH);
            
        end
        

        function [state, report] = updateState(model, state, problem, dx, drivingForces)

            [state, report] = updateState@AqueousModel(model, state, problem, dx, drivingForces);

            elyte = 'Electrolyte';
            
            % Cap quasi-particle concentrations of the "always positive" quasiparticle
            ics = model.(elyte).qpdict('Mg');           
            state.qpcs{ics} = max(0, state.qpcs{ics});
            ics = model.(elyte).qpdict('Cl');           
            state.qpcs{ics} = max(0, state.qpcs{ics});
            
        end
        
        
    end


end

