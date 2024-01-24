classdef PlatiniumCatalystLayer < CatalystLayer

    properties

    end
    
    methods
        
        function model = PlatiniumCatalystLayer(inputparams)

            model = model@CatalystLayer(inputparams);
            
        end
        
        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@CatalystLayer(model);

            % Assemble the reaction rate constants
            fn = @() PlatiniumCatalystLayer.updateReactionRateConstants;
            inputnames = {'cOHelyte', 'H2OaElyte', 'cOHinmr'};
            model = model.registerPropFunction({'elyteReactionRateConstant', fn, inputnames});
            model = model.registerPropFunction({'inmrReactionRateConstant', fn, inputnames});

            % Assemble equilibrium Potential for electrolyte
            fn = @() PlatiniumCatalystLayer.updateEelyte;
            inputnames = {'T', 'cOHelyte', 'pressureActiveGas', 'H2OaElyte'};
            model = model.registerPropFunction({'Eelyte', fn, inputnames});

            % Assemble equilibrium Potential for inmr
            fn = @() PlatiniumCatalystLayer.updateEinmr;
            inputnames = {'T', 'cOHinmr', 'pressureActiveGas', 'H2OaInmr'};
            model = model.registerPropFunction({'Einmr', fn, inputnames});

            % update source terms
            fn = @() PlatiniumCatalystLayer.updateSources;
            inputnames = {'elyteReactionRate', 'inmrReactionRate'};
            model = model.registerPropFunction({'elyteH2Osource', fn, inputnames});
            model = model.registerPropFunction({'inmrH2Osource', fn, inputnames});
            model = model.registerPropFunction({'elyteOHsource', fn, inputnames});
            model = model.registerPropFunction({'inmrOHsource', fn, inputnames});
            model = model.registerPropFunction({'activeGasSource', fn, inputnames});
            model = model.registerPropFunction({'eSource', fn, inputnames});

        end

        function state = updateReactionRateConstants(model, state)
            
            j0   = model.j0;
            cOH0 = model.sp.OH.c0;
            th   = 1e-2; % regularization parameter for the square root
            
            cOH = state.cOHelyte;
            aw  = state.H2OaElyte;
            
            state.elyteReactionRateConstant = j0.*regularizedSqrt(cOH/cOH0.*aw, th);
            
            aw = state.H2OaInmr;

            state.inmrReactionRateConstant = j0.*regularizedSqrt(aw, th);
            
        end
        
        function state = updateEelyte(model, state)

            T    = state.T;
            cOH  = state.cOHelyte;
            H2Oa = state.H2OaElyte;
            
            E0  = model.E0eff;
            c0  = model.sp.OH.c0;
            con = model.constants;

            F  = con.F;
            R  = con.R;

            state.Eelyte = E0 + R*T./(2*F).*log((H2Oa.^2).*(c0.^2).*(cOH.^-2));

        end
        
        function state = updateEinmr(model, state)

            T    = state.T;
            cOH  = state.cOHinmr;
            H2Oa = state.H2OaInmr;
            
            E0  = model.E0eff;
            c0  = model.sp.OH.c0;
            con = model.constants;

            F  = con.F;
            R  = con.R;

            state.Einmr = E0 + R*T./(2*F).*log((H2Oa.^2).*(c0.^2).*(cOH.^-2));
            
        end

        
        function state =  updateSources(model, state)

        % Reaction in Electrolyte : 2*H2O + 2*e- <<-> H2 + 2(OH-)_elyte
        % Reaction in Membrane :    2*H2O + 2*e- <<-> H2 + 2(OH-)_inmr
        % Here, the direction of the reaction that is indicated by the repeated arrow symbol corresponds to a positive computed reaction rate

            F = model.constants.F;
            vols = model.G.getVolumes();
            
            elyteR = state.elyteReactionRate;
            inmrR  = state.inmrReactionRate;

            R = elyteR + inmrR;
            
            state.activeGasSource = -R/F;
            state.elyteH2Osource  = 2*elyteR/F;
            state.inmrH2Osource   = 2*inmrR/F;
            state.elyteOHsource   = -2*elyteR/F;
            state.inmrOHsource    = -2*inmrR/F;
            state.eSource         = -2*R.*vols; % in the term eSource, we multiply by the volumes, see unit ([A])
            
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
