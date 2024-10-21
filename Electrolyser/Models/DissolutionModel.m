classdef DissolutionModel < BaseModel

    properties

        standardElectricalPotential                     % Standard electrical potential for the dissolution reaction [V]
        c0                     % Reference concentration [mol/m^3]
        referenceExchangeCurrentDensity                     % Reference exchange current density for the dissolution reaction [A/m^2]
        MW                     % Molar mass
        rho                    % Density
        volumeFraction0        % Initial volume fraction
        referenceVolumetricSurfaceArea % Initial volumetric surface area
        
        V         % Molar volume [m^3/mol]
        Np        % Number of particles [-]
        constants % physical constants
        
    end

    methods

        function model = DissolutionModel(inputparams)
            
            model = model@BaseModel();
            
            fdnames = {'G'              , ...
                       'standardElectricalPotential'             , ...
                       'c0'             , ...
                       'referenceExchangeCurrentDensity'             , ...
                       'MW'             , ...
                       'rho'            , ...
                       'volumeFraction0', ...
                       'referenceVolumetricSurfaceArea'};

            model = dispatchParams(model, inputparams, fdnames);

            model.constants = PhysicalConstants();

            % setup particule number
            vsa = model.referenceVolumetricSurfaceArea;
            vf  = model.volumeFraction0;
            model.Np = 1/(4*pi)*(vsa^3)/(3^2*vf^2);

            % setup molar volume
            model.V  = model.molecularWeight/model.rho;
            
        end

        function model = registerVarAndPropfuncNames(model)

        %% Declaration of the Dynamical Variables and Function of the model
        % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@BaseModel(model);

            varnames = {};
            % Temperature [K]
            varnames{end + 1} = 'T';
            % Volume Fraction [-]
            varnames{end + 1} = 'volumeFraction';
            % Reaction rate [A/m^3]
            varnames{end + 1} = 'reactionRate';
            % Electrical potential in electrolyte
            varnames{end + 1} = 'phiElyte';            
            % Electrical potential in solid
            varnames{end + 1} = 'phiSolid';
            % Equilibrium Electrical potential
            varnames{end + 1} = 'E';
            % Overpotential
            varnames{end + 1} = 'eta';
            % H+ concentration
            varnames{end + 1} = 'cH';
            % Volumetric surface area [m^2/m^3]
            varnames{end + 1} = 'volumetricSurfaceArea';
            % mass accumulation term (expressed in equivalent volume) [m^3/s]
            varnames{end + 1} = 'massAccum';
            % mass source term (expressed in equivalent volume) [m^3/s]
            varnames{end + 1} = 'massSource';
            % mass conservation equation (expressed in equivalent volume) [m^3/s]
            varnames{end + 1} = 'massCons';

            model = model.registerVarNames(varnames);

            fn = @() CatalystLayer.updateVolumetricSurfaceArea;
            inputnames = {'volumeFraction'};
            model = model.registerPropFunction({'volumetricSurfaceArea', fn, inputnames});

            fn = @() CatalystLayer.updateEta;
            inputnames = {'phiElyte', 'phiSolid', 'E'};
            model = model.registerPropFunction({'eta', fn, inputnames});
            
            fn = @() CatalystLayer.updateReactionRate;
            inputnames = {'volumetricSurfaceArea', 'eta', 'T'};
            model = model.registerPropFunction({'reactionRate', fn, inputnames});

            fn = @() CatalystLayer.updateE;
            inputnames = {'cH'};
            model = model.registerPropFunction({'E', fn, inputnames});
            
            fn = @() CatalystLayer.updateMassAccum;
            functionCallSetupFn = @(propfunction) PropFunction.accumFuncCallSetupFn(propfunction);
            fn = {fn, functionCallSetupFn};
            inputnames = {'volumeFraction'};
            model = model.registerPropFunction({'massAccum', fn, inputnames});

            fn = @() CatalystLayer.updateMassSource;
            inputnames = {'reactionRate'};
            model = model.registerPropFunction({'massSource', fn, inputnames});
            
            fn = @() CatalystLayer.updateMassCons;
            inputnames = {'massSource', 'massAccum'};
            model = model.registerPropFunction({'massCons', fn, inputnames});

        end


            function state = updateVolumetricSurfaceArea(model, state)

                Np = model.Np;

                vf = state.volumeFraction;

                a = (4*pi*Np)^(1/3)*(3^(2/3))*vf.^(2/3);

                state.volumetricSurfaceArea = a;
            end

            function state = updateEta(model, state)

                E        = state.E;
                phiElyte = state.phiElyte;                
                phiSolid = state.phiSolid;

                state.eta = phiSolid - phiElyte - E;
                
            end


            function state = updateReactionRate(model, state)

            % Reaction : IrO2 + 2 H2O <->> IrO4^2- + 4H+ + 2e-
            % (Here, the direction of the reaction that is indicated by the repeated arrow symbol corresponds to a positive computed reaction rate)
                
                j0 = model.referenceExchangeCurrentDensity;

                vsa = state.volumetricSurfaceArea;
                eta = state.eta;
                T   = state.T;

                state.reactionRate = vsa.*ButlerVolmerEquation(j0, 0.5, 2, eta, T);
                
            end


            function state = updateE(model, state)

                E0 = model.standardElectricalPotential;
                c0 = model.referenceConcentration;
                R  = model.constants.R;
                F  = model.constants.F;
                
                cH = state.cH;
                T  = state.T;
                
                state.E = E0 - R*T./(2*F).*log(((c0)./cH).^4);
                
            end


            function state = updateMassAccum(model, state, state0, dt)

                vols = model.G.getVolumes();
                
                vf  = state.volumeFraction;
                vf0 = state0.volumeFraction;

                state.massAccum = 1/dt*vols.*(vf - vf0);
                
            end


            function state = updateMassSource(model, state)

                V    = model.V;
                F    = model.constants.F;
                vols = model.G.getVolumes();
                
                rRate = state.reactionRate;

                state.massSource = -V/(2*F)*rRate.*vols;
                
            end

            function state = updateMassCons(model, state)

                source = state.massSource;
                accum  = state.massAccum;

                state.massCons = accum - source;
                
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
