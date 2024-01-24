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




%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
