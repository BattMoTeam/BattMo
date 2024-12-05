classdef HydrogenPorousTransportLayer < PorousTransportLayer

    properties

        % The species structure should now include
        % species.H2.molecularWeight : Molecular weight of H2 [kg/mol]
        
    end
    
    methods
        
        function model = HydrogenPorousTransportLayer(inputparams)

            model = model@PorousTransportLayer(inputparams);

            % add the H2 component in the indexing structures
            model.compInd.H2 = model.compInd.activeGas;
            model.gasInd.H2  = model.gasInd.activeGas;
            
            model.gasMW = model.species.H2.molecularWeight;
            
        end
        
        function model = registerVarAndPropfuncNames(model)
            
            model = registerVarAndPropfuncNames@PorousTransportLayer(model);

            phaseInd  = model.phaseInd;
            liquidInd = model.liquidInd;
            gasInd    = model.gasInd;
            compInd   = model.compInd;

            ncomp   = compInd.ncomp;
            ngas    = gasInd.ngas;
            nph     = phaseInd.nphase;
            nliquid = liquidInd.nliquid;
            nmobph  = numel(phaseInd.mobile);

            model = model.registerVarName('H2rhoeps');

            % assemble gas pressure using ideal gas law
            fn = @() HydrogenPorousTransportLayer.updateGasPressure;
            inputnames = {'H2rhoeps', 'H2Ogasrhoeps', 'T'};
            model = model.registerPropFunction({VarName({}, 'phasePressures', nph, phaseInd.gas), fn, inputnames});
            model = model.registerPropFunction({VarName({}, 'compGasPressures', ngas), fn, inputnames});

            % assemble H2mass
            fn = @() HydrogenPorousTransportLayer.updateH2GasMass;
            inputnames = {'H2rhoeps'};
            model = model.registerPropFunction({VarName({}, 'compGasMasses', ngas, gasInd.H2), fn, inputnames});
            
            fn = @() HydrogenPorousTransportLayer.updateGasViscosity;
            inputnames = {'T'};
            model = model.registerPropFunction({VarName({}, 'viscosities', nph, phaseInd.gas), fn, inputnames});
            
            
        end

        function state = updateH2GasMass(model, state)

            gasInd = model.gasInd;

            state.compGasMasses{gasInd.H2} = state.H2rhoeps;

        end
        
        function state = updateGasPressure(model, state)

            gasInd   = model.gasInd;
            phaseInd = model.phaseInd;
            MWH2O    = model.species.H2O.molecularWeight;
            MWH2     = model.species.H2.molecularWeight;
            R        = model.constants.R;
            
            mH2  = state.H2rhoeps;
            mH2O = state.H2Ogasrhoeps;
            vf   = state.volumeFractions{model.phaseInd.gas};
            T    = state.T;
            pH2  = R*T.*(mH2./MWH2)./vf;
            pH2O = R*T.*(mH2O./MWH2O)./vf;
            
            state.compGasPressures{gasInd.H2}  = pH2;
            state.compGasPressures{gasInd.H2O} = pH2O;
            state.phasePressures{phaseInd.gas} = pH2 + pH2O;
            
        end
        
        
        function state = updateGasViscosity(model, state)
            
            T = state.T;

            % Viscosity of H2 gas [Pa s]
            state.viscosities{model.mobPhaseInd.gas} = (1.644 + T.*(0.0278 - 1.17e-5.*T)).* 1e-6;

            
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
