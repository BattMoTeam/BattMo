classdef CatalystLayer < BaseModel

    properties

        constants

        compInd   % mapping structure for component indices
        phaseInd  % mapping structure for phase indices
        liquidInd % mapping structure for component indices
        gasInd    % mapping structure for component indices

        referenceExchangeCurrentDensity % Exchange current density
        standardElectricalPotential
        referencePotential
        E0eff
        species % species struct with field
        % - OH.chargeNumber  : Charge
        % - OH.referenceConcentration : OH reference concentration

        numberOfElectronsTransferred % Number of electron transfer
        
        chargeTransferCoefficient                 % coefficient in the exponent in Butler-Volmer equation [-]
        ionomerFractionArea                 % Fraction of specific area that is coversed with ionomer [-]
        volumetricSurfaceArea0 % Volumetric surface area [m^ -1]

        tortuosity % Tortuosity [-]

        include_dissolution
        DissolutionModel
        
    end

    methods

        function model = CatalystLayer(inputparams)
            
            model = model@BaseModel();
            
            fdnames = { 'G'                     , ...
                        'referenceExchangeCurrentDensity'                    , ...
                        'standardElectricalPotential'                    , ...
                        'referencePotential'                  , ...
                        'species'                    , ...
                        'numberOfElectronsTransferred'                     , ...
                        'chargeTransferCoefficient'                 , ...
                        'ionomerFractionArea'                 , ...
                        'volumetricSurfaceArea0', ...
                        'include_dissolution'   , ...
                        'tortuosity'};
            
            model = dispatchParams(model, inputparams, fdnames);

            if inputparams.include_dissolution
                model.DissolutionModel = DissolutionModel(inputparams.DissolutionModel);
            else
                model.subModelNameList = {};
            end
            
            model.E0eff = model.standardElectricalPotential - model.referencePotential;
            model.constants = PhysicalConstants();
            
        end

        function model = registerVarAndPropfuncNames(model)

        %% Declaration of the Dynamical Variables and Function of the model
        % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@BaseModel(model);

            varnames = {};
            % Temperature
            varnames{end + 1} = 'T';
            % CatalystLayer electrical potential (single value)
            varnames{end + 1} = 'E';
            % CatalystLayer total current (scalar value)
            varnames{end + 1} = 'I';
            % Electric Potential in electrolyte and ionomer
            varnames{end + 1} = 'phiElyte';
            varnames{end + 1} = 'phiInmr';
            % concentration of OH in electrotyte and ionomer
            varnames{end + 1} = 'cOHelyte';
            varnames{end + 1} = 'cOHinmr';
            % Equilibrium Potential for electrolyte and ionomer
            varnames{end + 1} = 'Eelyte';
            varnames{end + 1} = 'Einmr';
            % Partial pressure of the active gas (H2 or O2 for exampler) in the porous transport layer
            varnames{end + 1} = 'pressureActiveGas';
            % Potential difference with respect to the electrolyte and ionomer
            varnames{end + 1} = {'etaElyte'};
            varnames{end + 1} = {'etaInmr'};
            % Water activity in electrolyte and ionomer
            varnames{end + 1} = 'H2OaElyte';
            varnames{end + 1} = 'H2OaInmr';

            % Volumetric surface area [m^2/m^3]
            varnames{end + 1} = 'volumetricSurfaceArea';

            % Reaction rate constant (j0) for electrolyte and ionomer
            varnames{end + 1} = 'elyteReactionRateConstant';
            varnames{end + 1} = 'inmrReactionRateConstant';
            % Reaction rate for electrolyte and ionomer [A m^-3]
            varnames{end + 1} = 'elyteReactionRate';
            varnames{end + 1} = 'inmrReactionRate';
            
            % current source [A]. It is positive if there is a positive current source for the electrode: electrons are
            % removed in the catalyst layers (see redox equations). It is used to compute the control.
            varnames{end + 1} = 'eSource'; 
            % The following source terms are per volume. The units are [mol s^-1 m^-3]
            varnames{end + 1} = 'activeGasSource';
            varnames{end + 1} = 'elyteH2Osource';
            varnames{end + 1} = 'inmrH2Osource';
            varnames{end + 1} = 'elyteOHsource';
            varnames{end + 1} = 'inmrOHsource';
            
            model = model.registerVarNames(varnames);

            % Assemble equilibrium Potential for electrolyte
            fn = @() CatalystLayer.updateEelyte;
            inputnames = {'T', 'cOHelyte', 'pressureActiveGas', 'H2OaElyte'};
            model = model.registerPropFunction({'Eelyte', fn, inputnames});

            % Assemble equilibrium Potential for inmr
            fn = @() CatalystLayer.updateEinmr;
            inputnames = {'T', 'cOHinmr', 'pressureActiveGas', 'H2OaInmr'};
            model = model.registerPropFunction({'Einmr', fn, inputnames});

            % Assemble reactive potential
            fn = @() CatalystLayer.updateEtas;
            inputnames = {'phiElyte', 'phiInmr', 'Eelyte','Einmr', 'E'};
            model = model.registerPropFunction({'etaElyte', fn, inputnames});            
            model = model.registerPropFunction({'etaInmr', fn, inputnames});

            % Assemble the reaction rates
            fn = @() CatalystLayer.updateReactionRates;
            inputnames = {'volumetricSurfaceArea'    , ...
                          'elyteReactionRateConstant', ...
                          'etaElyte'                 , ...
                          'inmrReactionRateConstant' , ...
                          'etaInmr'};
            model = model.registerPropFunction({'elyteReactionRate', fn, inputnames});
            model = model.registerPropFunction({'inmrReactionRate', fn, inputnames});            

            fn = @() CatalystLayer.updateI;
            inputnames = {'eSource'};
            model = model.registerPropFunction({'I', fn, inputnames});

            if model.include_dissolution

                dm = 'DissolutionModel';
                
                fn = @() Catalyser.dispatchToDissolutionModel;
                inputnames = {'cOHelyte', 'phiElyte', 'E', 'T'};
                dissolutionVarnames = {'T', 'phiElyte', 'phiSolid', 'cH'};
                for idvar = 1 : numel(dissolutionVarnames)
                    varname = {dm, dissolutionVarnames{idvar}};
                    model = model.registerPropFunction({varname, fn, inputnames});
                end

                fn = @() Catalyser.updateVolumetricSurfaceAreaDissolution;
                inputnames = {{dm, 'volumetricSurfaceArea'}};
                model = model.registerPropFunction({'volumetricSurfaceArea', fn, inputnames});
                
            end

        end

        function state = updateVolumetricSurfaceAreaDissolution(model, state)

            dm = 'DissolutionModel';
            
            state.volumetricSurfaceArea = state.(dm).volumetricSurfaceArea;
            
        end
        
        
        function state = dispatchToDissolutionModel(model, state)

            dm = 'DissolutionModel';

            cOH      = state.cOHelyte;
            phiElyte = state.phiElyte;
            E        = state.E;
            T        = state.T;

            state.(dm).T = T;
            state.(dm).phiSolid = E;
            state.(dm).phiElyte = phiElyte;
            state.(dm).cH = 10^(-14)*((mol/litre)^2)./cOH;

        end

        
        function state = updateSources(model, state)
            error('virtual function. Should be overloaded by specific catalystlayer');
        end
        
        function state = updateEelyte(model, state)
            error('virtual function. Should be overloaded by specific catalystlayer');
        end

        function state = updateEinmr(model, state)
            error('virtual function. Should be overloaded by specific catalystlayer');            
        end

        function state = updateEtas(model, state)

            E        = state.E;
            phiElyte = state.phiElyte;
            phiInmr  = state.phiInmr;
            Einmr    = state.Einmr;
            Eelyte   = state.Eelyte;

            state.etaElyte = E - phiElyte - Eelyte;
            state.etaInmr  = E - phiInmr - Einmr;
            
        end

        function state = updateReactionRateConstants(model, state)

            error('implemented in derived class')
            
        end
            
        function state = updateReactionRates(model, state)

            Xinmr = model.ionomerFractionArea;
            alpha = model.chargeTransferCoefficient;
            n     = model.numberOfElectronsTransferred;

            vsa      = state.volumetricSurfaceArea;
            etaElyte = state.etaElyte;
            etaInmr  = state.etaInmr;
            j0elyte  = state.elyteReactionRateConstant;
            j0inmr   = state.inmrReactionRateConstant;
            T        = state.T;
            
            jElyte = (1 - Xinmr)*vsa.*ButlerVolmerEquation(j0elyte, alpha, n, etaElyte, T);
            jInmr  = Xinmr*vsa.*ButlerVolmerEquation(j0inmr, alpha, n, etaInmr, T);

            state.elyteReactionRate = jElyte;
            state.inmrReactionRate  = jInmr;
            
        end

        function state = updateI(model, state)

            state.I = sum(state.eSource);
            
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
