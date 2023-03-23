classdef PorousTransportLayerBoundary < BaseModel
    
    properties

        compInd     % mapping structure for component indices
        phaseInd    % mapping structure for phase indices
        mobPhaseInd % mapping structure for mobile phase indices
        liquidInd   % mapping structure for component indices
        gasInd      % mapping structure for component indices

        controlValues % Structure with following fields
                      % - cOH
                      % - liqrho
                      % - mobilePhasePressures (only for mobile phases)
                      % - gasDensities
        
    end
    
    methods
        
        function model = PorousTransportLayerBoundary(paramobj)
            
            model = model@BaseModel();

            fdnames = {'compInd'    , ...
                       'phaseInd'   , ...
                       'liquidInd'  , ...
                       'mobPhaseInd', ...
                       'gasInd'};
            model = dispatchParams(model, paramobj, fdnames);

            
        end

        function model = registerVarAndPropfuncNames(model)
            
        %% Declaration of the Dynamical Variables and Function of the model
        % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@BaseModel(model);

            phaseInd  = model.phaseInd;
            liquidInd = model.liquidInd;
            gasInd    = model.gasInd;
            compInd   = model.compInd;

            ncomp   = compInd.ncomp;
            ngas    = gasInd.ngas;
            nph     = phaseInd.nphase;
            nliquid = liquidInd.nliquid;
            nmobph  = numel(phaseInd.mobile);

            varnames = {};

            % Phase pressures in [Pa]
            phasePressures = VarName({}, 'phasePressures', nph);
            varnames{end + 1} = phasePressures;

            % Gas densities in mass per gas volume (one value for each gas component)
            gasDensities = VarName({}, 'gasDensities', ngas);
            varnames{end + 1} = gasDensities;

            % Liquid density (Mass of liquid per volume of liquid) in [kg m^-3]
            varnames{end + 1} = 'liqrho';

            varnames{end + 1} = 'cOH';

            % electric potential
            % not needed because (for the moment) we impose a no current condition at the boundary
            % varnames{end + 1} = 'phi';

            % Boundary equations
            bcEquations = VarName({}, 'bcEquations', 2 + ngas);
            varnames{end + 1} = bcEquations;
            
            % Boundary control equations
            bcControlEquations = VarName({}, 'bcControlEquations', nmobph);
            varnames{end + 1} = bcControlEquations;

            model = model.registerVarNames(varnames);

            model = model.removeVarName(VarName({}, 'phasePressures', nph, phaseInd.solid));
            
            fn = @() PorousTransportLayerBoundary.setupBcControlEquations;
            inputnames = {VarName({}, 'phasePressures', nph, model.mobPhaseInd.phaseMap)};
            model = model.registerPropFunction({bcControlEquations, fn, inputnames});
            
        end

        function state = setupBcControlEquations(model, state)

            nmobph = model.mobPhaseInd.nmobphase;
            
            for imobph = 1 : nmobph
                iph = model.mobPhaseInd.phaseMap(imobph);
                eqs{imobph} = state.phasePressures{iph} - model.controlValues.mobilePhasePressures{imobph};
            end

            state.bcControlEquations = eqs;
            
        end
        
        
    end


end




%{
Copyright 2021-2023 SINTEF Industry, Sustainable Energy Technology
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
