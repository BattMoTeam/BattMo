classdef SEIActiveMaterial < ActiveMaterial

    properties

        SolidElectrodeInterface
        SideReaction

    end

    methods

        function model = SEIActiveMaterial(inputparams)

            model = model@ActiveMaterial(inputparams);

            model.SideReaction = SideReaction(inputparams.SideReaction);
            model.SolidElectrodeInterface = SolidElectrodeInterface(inputparams.SolidElectrodeInterface);

        end

        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@ActiveMaterial(model);

            itf = 'Interface';
            sd  = 'SolidDiffusion';
            sei = 'SolidElectrodeInterface';
            sr  = 'SideReaction';

            varnames = {};
            % Total reaction rate
            varnames{end + 1} = 'totalReactionRate';
            % SEI charge consservation equation
            varnames{end + 1} = 'seiInterfaceChargeCons';
            % SEI charge consservation equation
            varnames{end + 1} = {itf, 'externalPotentialDrop'};

            model = model.registerVarNames(varnames);

            model = model.setAsStaticVarName({sei, 'cExternal'});
            
            if model.isRootSimulationModel
                model = model.setAsStaticVarName({sr, 'phiElectrolyte'});
            end

            fn = @SEIActiveMaterial.assembleSEIchargeCons;
            inputnames = {{itf, 'intercalationFlux'}, {sr, 'reactionRate'}, 'totalReactionRate'};
            model = model.registerPropFunction({'seiInterfaceChargeCons', fn, inputnames});

            fn = @SEIActiveMaterial.dispatchTemperature;
            model = model.registerPropFunction({{sr, 'T'}, fn, {'T'}});

            fn = @SEIActiveMaterial.updatePotentialDrop;
            inputnames = {{'totalReactionRate'}, {'SolidElectrodeInterface', 'delta'}};
            model = model.registerPropFunction({{itf, 'externalPotentialDrop'}, fn, inputnames});
            model = model.registerPropFunction({{sr, 'externalPotentialDrop'}, fn, inputnames});

            fn = @SEIActiveMaterial.Interface.updateEtaWithEx;
            % Comment about the syntax used below : Here we use the more cumbersome syntax (using VarName)
            % because we use a different model than the current model in the definition of the propfunction. The
            % handy-syntax could have been implemented here too (but this has not been done...)
            varname = VarName({itf}, 'eta');
            inputvarnames = {VarName({itf}, 'phiElectrolyte'), ...
                             VarName({itf}, 'phiElectrode')  , ...
                             VarName({itf}, 'OCP')           , ...
                             VarName({itf}, 'externalPotentialDrop')};
            modelnamespace = {itf};
            propfunction = PropFunction(varname, fn, inputvarnames, modelnamespace);
            model = model.registerPropFunction(propfunction);

            fn = @SEIActiveMaterial.dispatchSEIRate;
            inputnames = {{sr, 'reactionRate'}};
            model = model.registerPropFunction({{sei, 'reactionRate'}, fn, inputnames});

            fn = @SEIActiveMaterial.updateSEISurfaceConcentration;
            inputnames = {{sei, 'cInterface'}};
            model = model.registerPropFunction({{sr, 'c'}, fn, inputnames});

            if model.isRootSimulationModel

                fn = @SEIActiveMaterial.updatePhi;
                model = model.registerPropFunction({{sr, 'phiElectrode'}, fn, {'E'}});

                fn = @SEIActiveMaterial.updateChargeCons;
                inputnames = {'I', ...
                              'totalReactionRate'};
                model = model.registerPropFunction({'chargeCons', fn, inputnames});

            end

        end


        function model = setupForSimulation(model)

            model = model.equipModelForComputation();

            itf = 'Interface';
            sd  = 'SolidDiffusion';
            sei = 'SolidElectrodeInterface';
            sr  = 'SideReaction';

            % add scaling 
            F = model.(sd).constants.F;
            scalingcoef = F;

            scalings = {{{'chargeCons'}, scalingcoef}};

            % add scaling
            rp  = model.(sd).particleRadius;
            vsa = model.(sd).volumetricSurfaceArea;
            scalingcoef = vsa*(4*pi*rp^3/3);

            scalings = horzcat(scalings, ...
                               {{{sd, 'massCons'}, scalingcoef}, ...
                                {{sd, 'solidDiffusionEq'}, scalingcoef}});

            % add scaling
            
            Mw  = model.(sei).molecularWeight;
            rho = model.(sei).density;
            scalingcoef = 0.5*Mw/rho;

            scalings = horzcat(scalings, ...
                               {{{sei, 'widthEq'}, scalingcoef}});

            % add scaling

            deltaref = 1*nano*meter;
            scalingcoef = deltaref;
            
            scalings = horzcat(scalings, ...
                               {{{sei, 'interfaceBoundaryEq'}, scalingcoef}, ...
                                {{sei, 'massCons'}, scalingcoef}});


            model.scalings = scalings;
            
        end

        function state = dispatchTemperature(model, state)

            state = dispatchTemperature@ActiveMaterial(model, state);
            state.SideReaction.T = state.T;

        end

        function state = dispatchE(model, state)

            state.phi = state.E;

        end

        function state = setupControl(model, state)

            state.totalReactionRate = state.controlCurrentSource;

        end

        function state = updateSEISurfaceConcentration(model, state)

            sei = 'SolidElectrodeInterface';
            sr  = 'SideReaction';

            state.(sr).c = state.(sei).cInterface;

        end

        function state = updateChargeCons(model, state)
            itf = 'Interface';
            F    = model.(itf).constants.F;

            I = state.I;
            R = state.totalReactionRate;

            state.chargeCons = I - R*F;

        end

        function state = dispatchSEIRate(model, state)

            sei = 'SolidElectrodeInterface';
            sr  = 'SideReaction';

            state.(sei).reactionRate = state.(sr).reactionRate;

        end

        function state = updatePhi(model, state)

            state = updatePhi@ActiveMaterial(model, state);
            state.SideReaction.phiElectrode = state.E;

        end

        function state = assembleSEIchargeCons(model, state)

            itf = 'Interface';
            sr  = 'SideReaction';

            vsa = model.(itf).volumetricSurfaceArea;

            Rint = state.(itf).intercalationFlux;
            Rsei = state.(sr).reactionRate;
            R    = state.totalReactionRate;

            state.seiInterfaceChargeCons = R - Rint - Rsei;

        end

        function state = updatePotentialDrop(model, state)

            sei = 'SolidElectrodeInterface';
            itf = 'Interface';
            sr  = 'SideReaction';

            kappaSei = model.(sei).conductivity;
            F = model.(itf).constants.F;

            R     = state.totalReactionRate;
            delta = state.(sei).delta;

            dphi = F*delta.*R/kappaSei;

            state.(itf).externalPotentialDrop = dphi;
            state.(sr).externalPotentialDrop  = dphi;

        end

        function cleanState = addStaticVariables(model, cleanState, state, state0)

            cleanState = addStaticVariables@ActiveMaterial(model, cleanState, state);

            sei = 'SolidElectrodeInterface';
            sr  = 'SideReaction';

            cleanState.(sei).cExternal     = state.(sei).cExternal;

            if model.isRootSimulationModel
                
                cleanState.(sr).phiElectrolyte = state.(sr).phiElectrolyte;

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
