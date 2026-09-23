function [energy, extras] = advancedComputeCellEnergy(model, varargin)

    opt = struct('verbose', false);
    opt = merge_options(opt, varargin{:});
    
    if isstruct(model)
        model = setupModelFromJson(model);
    end

    ecs = EquilibriumConcentrationSolver(model);

    % compute initial guess
    [stateInit, ecs] = ecs.setupInitialState();
    stateInit = ecs.evalVarName(stateInit, 'voltage');
    stateInit = ecs.computeConcentrations(stateInit.voltage, 'verbose', opt.verbose);;

    [stateStart, stateEnd] = ecs.computeExtremalStates(stateInit, 'verbose', opt.verbose);
    
    N = 100;

    theta = linspace(0, 1, N);

    ne = 'NegativeElectrode';
    pe = 'PositiveElectrode';

    eldes = {ne, pe};

    for ielde = 1 : numel(eldes)

        elde = eldes{ielde};
        
        nam = ecs.(elde).numberOfActiveMaterial;

        theta_elde = {};
        
        for iam = 1 : nam

            theta_elde{iam} = linspace(stateStart.(elde).stoichiometries{iam}, ...
                                       stateEnd.(elde).stoichiometries{iam}  , ...
                                       N);
            state.(elde).stoichiometries{iam} = theta_elde{iam};
            
        end

        thetas.(elde) = theta_elde;
        
    end

    state = ecs.evalVarName(state, 'ocps');

    F = PhysicalConstants.F;

    E = 0;
    
    for ielde = 1 : numel(eldes)
        elde = eldes{ielde};
        nam = ecs.(elde).numberOfActiveMaterial;
        for iam = 1 : nam
            
            theta = thetas.(elde){iam};
            ocp = state.(elde).ocps{iam};
            ocp = 0.5*(ocp(1 : end - 1) + ocp(2 : end));
            dtheta = diff(theta);

            E = E + F*ecs.(elde).volumes(iam)*sum(ocp.*dtheta)*ecs.(elde).saturationConcentrations(iam);
            
        end
        
    end

    energy = E;

    if nargout > 1

        stateStart = ecs.evalVarName(stateStart, 'voltage');
        Ustart = stateStart.voltage;
        
        stateEnd = ecs.evalVarName(stateEnd, 'voltage');
        Uend = stateEnd.voltage;

        N = 50;
        
        U = linspace(Ustart, Uend, N)';

        [state, failure, ecs] = ecs.computeConcentrations(U, 'verbose', opt.verbose);
        
        if failure
            error('could not find solution');
        end
        
        state = ecs.evalVarName(state, {pe, 'amount'});

        a = state.(pe).amount;
        soc = (a - a(1))./(a(end) - a(1));

        dischargeFunction = @(s) interp1(soc, U, s, 'linear');

        extras.dischargeFunction = dischargeFunction;
        
    end
    
    
end

%{
Copyright 2021-2026 SINTEF Industry, Sustainable Energy Technology
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
