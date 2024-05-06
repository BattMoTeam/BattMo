classdef Interface < BaseModel

    properties

        % Physical constants
        constants = PhysicalConstants();


        %% Input parameters

        % Standard parameters

        saturationConcentration      % the saturation concentration of the guest molecule in the host material
        numberOfElectronsTransferred % stoichiometric number of electrons transferred in the electrochemical reaction
        volumetricSurfaceArea        % surface area of the active material - electrolyte interface per volume of electrode
        activationEnergyOfReaction   % the activation energy of the electrochemical reaction
        reactionRateConstant         % the reaction rate constant of the electrochemical reaction

        % The exchange current at the electrode-electrolyte interface under equilibrium conditions.
        % Tf empty, the default expression using the reaction rate constant is used, see method
        % Interface.updateReactionRateCoefficient. The function is given as a struct with the fields:
        %   - type = {"function", "constant"} % if "constant" is selected, we use the reactionRateConstant value
        %   - functionname :  matlab function name (should be available in path)
        %   - argumentlist = ["cElectrodeSurface", "cmax"]
        exchangeCurrentDensity

        guestStoichiometry100 % the ratio of the concentration of the guest molecule to the saturation concentration
                              % of the guest molecule in a phase at a cell voltage that is defined as 100% SOC
        guestStoichiometry0   % the ratio of the concentration of the guest molecule to the saturation concentration
                              % of the guest molecule in a phase at a cell voltage that is defined as 0% SOC
        density               % the mass density of the active material

        % A function to determine the open-circuit potential of the electrode under given conditions
        %   - type : "function";
        %   - functionname :  matlab function name (should be available in path)
        %   - argumentlist : ["cElectrode", "T", "cmax"]
        openCircuitPotential

        chargeTransferCoefficient % the charge transfer coefficient that enters in the Butler-Volmer equation (symbol: alpha)

        sei_type % string defining interface type. Can take value
                 % - 'none' (default)
                 % - 'bolay'

        %% Computed parameters at model setup

        computeOCPFunc % Function handler to compute OCP
        useJ0Func      % true if we use a function to compute the function computeJ0Func to compute the exchange current density
        computeJ0Func  % used when useJ0Func is true. Function handler to compute J0 as function of cElectrode, see method updateReactionRateCoefficient

    end

    methods

        function model = Interface(inputparams)

            model = model@BaseModel();

            fdnames = {'G'                           , ...
                       'saturationConcentration'     , ...
                       'numberOfElectronsTransferred', ...
                       'volumetricSurfaceArea'       , ...
                       'activationEnergyOfReaction'  , ...
                       'reactionRateConstant'        , ...
                       'exchangeCurrentDensity'      , ...
                       'guestStoichiometry100'       , ...
                       'guestStoichiometry0'         , ...
                       'density'                     , ...
                       'openCircuitPotential'        , ...
                       'chargeTransferCoefficient'   , ...
                       'sei_type'};

            model = dispatchParams(model, inputparams, fdnames);

            model.computeOCPFunc = str2func(inputparams.openCircuitPotential.functionname);

            j0 = inputparams.exchangeCurrentDensity;

            if ~isempty(j0)
                switch j0.type
                  case 'function'
                    model.useJ0Func = true;
                    model.computeJ0Func = str2func(j0.functionname);
                  case 'constant'
                    model.useJ0Func = false;
                  otherwise
                    error('type of j0 not recognized.')
                end
            else
                model.useJ0Func = false;
            end

        end

        function model = registerVarAndPropfuncNames(model)

            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@BaseModel(model);

            varnames = {};
            % Temperature
            varnames{end + 1} = 'T';
            % potential in electrode
            varnames{end + 1} = 'phiElectrode';
            % charge carrier concentration in electrode - value at surface
            varnames{end + 1} = 'cElectrodeSurface';
            % potential in electrolyte
            varnames{end + 1} = 'phiElectrolyte';
            % charge carrier concentration in electrolyte
            varnames{end + 1} = 'cElectrolyte';
            % Electrode over potential
            varnames{end + 1} = 'eta';
            % Reaction rate [mol s^-1 m^-2]
            varnames{end + 1} = 'R';
            % External potential drop used in Butler-Volmer in case of SEI see :class:`Electrochemistry.SEIActiveMaterial`
            % varnames{end + 1} = 'externalPotentialDrop';
            %
            varnames{end + 1} = 'dUdT';
            % OCP [V]
            varnames{end + 1} = 'OCP';
            % Reaction rate coefficient [A m^-2]
            varnames{end + 1} = 'j0';

            model = model.registerVarNames(varnames);

            fn = @Interface.updateReactionRateCoefficient;
            if model.useJ0Func
                inputnames = {'cElectrodeSurface'};
            else
                inputnames = {'T', 'cElectrolyte', 'cElectrodeSurface'};
            end
            model = model.registerPropFunction({'j0', fn, inputnames});

            fn = @Interface.updateOCP;
            inputnames = {'cElectrodeSurface', 'T'};
            model = model.registerPropFunction({'OCP', fn, inputnames});
            model = model.registerPropFunction({'dUdT', fn, inputnames});

            fn = @Interface.updateEta;
            inputnames = {'phiElectrolyte', 'phiElectrode', 'OCP'};
            model = model.registerPropFunction({'eta', fn, inputnames});

            % This function is used when SEI layer
            % fn = @Interface.updateEtaWithEx;
            % inputnames = {'phiElectrolyte', 'phiElectrode', 'OCP', 'externalPotentialDrop'};
            % model = model.registerPropFunction({'eta', fn, inputnames});

            fn = @Interface.updateReactionRate;
            inputnames = {'T', 'eta', 'j0'};
            model = model.registerPropFunction({'R', fn, inputnames});


        end

        function jsonstruct = exportParams(model)

            jsonstruct = exportParams@BaseModel(model);

            fdnames = {'saturationConcentration'     , ...
                       'numberOfElectronsTransferred', ...
                       'volumetricSurfaceArea'       , ...
                       'activationEnergyOfReaction'  , ...
                       'reactionRateConstant'        , ...
                       'exchangeCurrentDensity'      , ...
                       'guestStoichiometry100'       , ...
                       'guestStoichiometry0'         , ...
                       'openCircuitPotential'        , ...
                       'chargeTransferCoefficient'};

            for ifd = 1 : numel(fdnames)
                fdname = fdnames{ifd};
                jsonstruct.(fdname) = model.(fdname);
            end


        end

        function state = dipatchTemperature(model, state)

            sd = 'SolidDiffusion';
            state.(sd).T = state.T;

        end

        function state = updateOCP(model, state)

            computeOCP = model.computeOCPFunc;
            cmax = model.saturationConcentration;

            c = state.cElectrodeSurface;
            T = state.T;

            [state.OCP, state.dUdT] = computeOCP(c, T, cmax);

        end

        function state = updateReactionRateCoefficient(model, state)


            if model.useJ0Func

                computeJ0 = model.computeJ0Func;
                cmax      = model.saturationConcentration;
                theta0    = model.guestStoichiometry0;
                theta100  = model.guestStoichiometry100;

                c = state.cElectrodeSurface;

                cmin = theta0*cmax;
                cmax = theta100*cmax;

                soc = (c - cmin)./(cmax - cmin);

                j0 = computeJ0(soc);

            else

                Tref = 298.15;  % [K]

                cmax = model.saturationConcentration;
                k0   = model.reactionRateConstant;
                Eak  = model.activationEnergyOfReaction;
                n    = model.numberOfElectronsTransferred;
                F    = model.constants.F;
                R    = model.constants.R;

                T      = state.T;
                cElyte = state.cElectrolyte;
                c      = state.cElectrodeSurface;

                % Calculate reaction rate constant
                k = k0.*exp(-Eak./R.*(1./T - 1/Tref));

                % We use regularizedSqrt to regularize the square root function and avoid the blow-up of derivative at zero.
                th = 1e-3*cmax;
                coef = cElyte.*(cmax - c).*c;
                coef(coef < 0) = 0;
                j0 = k.*regularizedSqrt(coef, th)*n*F;

            end

            state.j0 = j0;

        end

        function state = updateEta(model, state)

            phiElyte = state.phiElectrolyte;
            phiElde  = state.phiElectrode;
            OCP      = state.OCP;

            state.eta = (phiElde - phiElyte - OCP);

        end

        function state = updateEtaWithEx(model, state)

            phiElyte = state.phiElectrolyte;
            phiElde  = state.phiElectrode;
            OCP      = state.OCP;
            dphi     = state.externalPotentialDrop;

            state.eta = (phiElde - phiElyte - OCP - dphi);

        end


        function state = updateReactionRate(model, state)
        % From definition of the overpotential eta, we have that reaction rate R is positive for oxydation.

            n     = model.numberOfElectronsTransferred;
            F     = model.constants.F;
            alpha = model.chargeTransferCoefficient;

            T   = state.T;
            j0  = state.j0;
            eta = state.eta;

            R = ButlerVolmerEquation(j0, alpha, n, eta, T);

            state.R = R/(n*F); % reaction rate in mol/(s*m^2)

        end

    end
end

%% References
%   [1] Torchio et al, Journal of The Electrochemical Society, 163 (7)
%   A1192-A1205 (2016), DOI: 10.1149/2.0291607jes



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
