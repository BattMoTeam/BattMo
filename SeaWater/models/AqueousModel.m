classdef AqueousModel < BaseModel

    properties

        Electrolyte
        totals
        totalnames
        pH
        
    end
    
    methods
        
        function model = AqueousModel(elyteModel, totals, totalnames, pH)

            model = model@BaseModel();
            
            model.nonlinearTolerance = 1e-16;
            
            model.Electrolyte = elyteModel;
            model.totals      = totals;
            model.totalnames  = totalnames;
            model.pH          = pH;
            
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            
            opts = struct('ResOnly', false, 'iteration', 0); 
            opts = merge_options(opts, varargin{:});

            elyte    = model.Electrolyte;
            tots     = model.totals;
            totnames = model.totalnames;
            pH       = model.pH;
            
            logspdict  = elyte.logspdict;
            qpdict     = elyte.qpdict;
            indsolidsp = elyte.indsolidsp;
            nsolid     = elyte.nsolid;
            
            % initialize AD
            if (not(opts.ResOnly))
                state = model.initStateAD(state);
            end

            state = elyte.updateFromLogConcentration(state);
            state = elyte.updateConcentrations(state);
            
            nc = size(state.qpcs{1}, 1);
            for isolid = 1 : nsolid
                indsp = indsolidsp(isolid);
                state.cs{indsp} = zeros(nc, 1);
            end
            
            % use Electrolyte model with volumeFraction = 1
            state.volumeFraction = 1;
            for iqp = 1 : elyte.nqp
                state.qpepscs{iqp} = state.qpcs{iqp};
            end
            state = elyte.updateAtomicMassCons(state);
            
            eqs = state.atomicMassCons;
            for ind = 1 : numel(totnames)
                nm = totnames{ind};
                indqp = qpdict(nm);
                eqs{end + 1} = state.qpcs{indqp} - tots(:, ind);
            end
            eqs{end + 1} = state.pcs{logspdict('H+')} - log(10^(-pH)*mol/litre);
            
            primaryVars = model.getPrimaryVariableNames();
            
            nc = numel(eqs);
            types = repmat({'cells'}, 1, nc);
            names = arrayfun(@(ind) sprintf('eqs_%d', ind), (1 : nc), 'uniformoutput', false);
        
            %% setup LinearizedProblem that can be processed by MRST Newton API
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
            
        end
        
        function  p = getPrimaryVariableNames(model, state)

            qp = arrayfun(@(i) {'qpcs', i}, (1 : model.Electrolyte.nqp)', 'uniformoutput', false);
            pcs = arrayfun(@(i) {'pcs', i}, (1 : model.Electrolyte.nlogsp)', 'uniformoutput', false);
            
            p = {qp{:}, pcs{:}};
            
        end       

        function [state, report] = updateState(model, state, problem, dx, drivingForces)

            [state, report] = updateState@BaseModel(model, state, problem, dx, drivingForces);

            elyte = 'Electrolyte';
            
            % Cap quasi-particle concentrations of the "always positive" quasiparticle
            % ics = model.(elyte).qpdict('Mg');           
            % state.qpcs{ics} = max(0, state.qpcs{ics});
            % ics = model.(elyte).qpdict('Cl');           
            % state.qpcs{ics} = max(0, state.qpcs{ics});
            
        end
        

        function state = solveAqueousMixture(model, stateguess)

            solver = NonLinearSolver();

            %% launch solver
            dt            = 0; 
            drivingForces = []; 

            [state, failure, report] = solveMinistep(solver, model, stateguess, stateguess, dt, drivingForces);
            
            converged = report.Converged;
            if ~converged
                nlreport = report.NonlinearReport{end};
                failure  = nlreport.Failure;
                if failure
                    msg = ['Model step resulted in failure state. Reason: ', ...
                           nlreport.FailureMsg]; %#ok<AGROW>
                    error(msg);
                end
            end
            
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
