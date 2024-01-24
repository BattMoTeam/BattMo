classdef MagnesiumBatteryNoPrecipitation < MagnesiumBattery

    methods

        function model = MagnesiumBatteryNoPrecipitation(inputparams)

            model = model@MagnesiumBattery(inputparams);

        end
        
        function model = registerVarAndPropfuncNames(model)
                    
            model = registerVarAndPropfuncNames@MagnesiumBattery(model);
            
        end

        function [problem, state] = getEquations(model, state0, state,dt, drivingForces, varargin)
                    
            opts = struct('ResOnly', false, 'iteration', 0);
            opts = merge_options(opts, varargin{:});

            time = state0.time + dt;
            if(not(opts.ResOnly))
                state = model.initStateAD(state);
            end

            elyte = 'Electrolyte';
            anam  = 'AnodeActiveMaterial';
            ctam  = 'CathodeActiveMaterial';
            ct    = 'Cathode';
            an    = 'Anode';

            % Initialize temperature (given in model)
            state = model.initializeTemperature(state);

            % Dispatch the temperature in all the submodels.
            state = model.dispatchTemperature(state);

            % Update the concentration of the main species from their logarithm values
            state.(elyte) = model.(elyte).updateFromLogConcentration(state.(elyte));

            % Update the concentration of the main species from their logarithm values
            state.(elyte) = model.(elyte).updateSolidVolumeFraction(state.(elyte));
            
            % update potential in anode
            state.(an) = model.(an).updatePotential(state.(an));

            % Dispatch concentration from electrolyte to anode active material
            state = model.updateAnodeActiveMaterialConcentration(state);

            % Update ENernst coefficient (anode active material)
            state.(anam) = model.(anam).updateENernst(state.(anam));

            % Update electrical potentials in active material from electrolyte and anode 
            % The updateAnodeActiveMaterialPotentials takes also care of updating the ASurfElectrode property in the anode
            state = model.updateAnodeActiveMaterialPotentials(state);

            % Update eta value in anode active material
            state.(anam) = model.(anam).updateEta(state.(anam));

            % Update Reaction rate in anode active material
            state.(anam) = model.(anam).updateReactionRate(state.(anam));

            % Update source term in anode
            state = model.updateAnodeSource(state);

            % Update accumulation term in anode
            state = model.updateAnodeAccumTerm(state, state0, dt);

            % Assemble mass conservation equation in anode
            state.(an) = model.(an).updateMassCons(state.(an));

            % update potential in cathode (we use infinite conductivity assumption for the moment)
            state.(ct) = model.(ct).updatePotential(state.(ct));

            % Dispatch concentration from electrolyte to cathode active material
            state = model.updateCathodeActiveMaterialConcentration(state);

            % Update ENernst coefficient (anode active material)
            state.(ctam) = model.(ctam).updateENernst(state.(ctam));

            % Update electrical potentials in cathode active material from electrolyte and cathode
            state = model.updateCathodeActiveMaterialPotentials(state);

            % Update eta value in cathode active material
            state.(ctam) = model.(ctam).updateEta(state.(ctam));

            % Update Reaction rate in cathode active material
            state.(ctam) = model.(ctam).updateReactionRate(state.(ctam));

            % Update source term in cathode
            state = model.updateCathodeSource(state);

            % Update accumulation term in cathode
            state = model.updateCathodeAccumTerm(state, state0, dt);

            % Assemble mass conservation equation in cathode
            state.(ct) = model.(ct).updateMassCons(state.(ct));

            % Update electrolyte current source (using reaction rates)
            state = model.updateElectrolyteCurrentSource(state);

            % Update electrolyte volume fraction (should some to one with the other volume fractions)
            state = model.updateElectrolyteVolumeFractionEquation(state);

            % Update the concentrations of the solutes (using chemistry)
            state.(elyte) = model.(elyte).updateConcentrations(state.(elyte));
            
            % Update the transference numbers from the concentration values
            state.(elyte) = model.(elyte).updateTransferenceNumbers(state.(elyte));

            % Update the current values
            state.(elyte) = model.(elyte).updateCurrent(state.(elyte));

            % Setup the charge conservation equation
            state.(elyte) = model.(elyte).updateChargeConservation(state.(elyte));

            % Update the diffusion fluxes for the quasi-particles
            state.(elyte) = model.(elyte).updateDiffFluxes(state.(elyte));

            % Update the migration fluxes for the quasi-particles
            state.(elyte) = model.(elyte).updateMigFluxes(state.(elyte));

            % Update the total fluxes for the quasi-particles (sum of the diffustion and migration fluxes)
            state.(elyte) = model.(elyte).updateFluxes(state.(elyte));

            % Update the quasi-particle source term (contributions from reactions at the electrodes and the solid
            % precipitation)
            state = model.updateQuasiParticleSource(state);

            % Update the electrolyte accumulation terms for the quasi-particles
            state = model.updateElectrolyteAccumTerms(state, state0, dt);

            % Assemble the mass conservation equations for the quasi-particles
            state.(elyte) = model.(elyte).updateMassCons(state.(elyte));

            % Assemble the atomic mass conservation equations (basically the definitions of the quasi-particles)
            state.(elyte) = model.(elyte).updateAtomicMassCons(state.(elyte));

            %% galvanostatic equation

            I = drivingForces.src(time);
            R = state.(ctam).R;
            F = model.con.F;
            vols = model.(ct).G.getVolumes();
            ctGalvanostatic = I + sum(2*F*R.*vols);

            % at the anode we impose zero potential
            anGalvanostatic = state.(an).E;
            
            %% We collect mass and charge conservation equations for the electrolyte and the electrodes
            nqp = model.(elyte).nqp;
            
                        
            anvols    = model.(an).G.getVolumes();
            ctvols    = model.(an).G.getVolumes();
            elytevols = model.(elyte).G.getVolumes();
            F = model.con.F;
            
            eqs = {};
            eqs{end + 1} = dt./anvols.*state.(an).massCons;
            eqs{end + 1} = dt./ctvols.*state.(ct).massCons;
            for ind = 1 : nqp
                eqs{end + 1} = dt./elytevols.*state.(elyte).qpMassCons{ind};
            end
            for ind = 1 : nqp
                eqs{end + 1} = 1e-3*state.(elyte).atomicMassCons{ind};
            end
            eqs{end + 1} = 1/F*state.(elyte).chargeCons;
            eqs{end + 1} = ctGalvanostatic;
            eqs{end + 1} = anGalvanostatic;
            eqs{end + 1} = state.(elyte).volumeFractionEquation;
            
            %% setup equation names (just for printed output)
            names = {'anode_masscons', 'cathode_masscons'};
            for ind = 1 : nqp
                names{end + 1} = sprintf('qp_%d_masscons', ind);
            end
            for ind = 1 : nqp
                names{end + 1} = sprintf('qp_%d_atommasscons', ind);
            end
            names = {names{:}, 'elyte_chargecons', 'ct_galvanostatic', 'an_galvanostatic', 'electrolyte_vf'};
            types = repmat({'cell'}, 1, numel(names));

            primaryVars = model.getPrimaryVariableNames();

            %% setup LinearizedProblem that can be processed by MRST Newton API
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

        end

        function p = getPrimaryVariableNames(model)
            
            qp = arrayfun(@(i) {'Electrolyte', 'qpepscs', i}, (1 : model.Electrolyte.nqp)', 'uniformoutput', false);
            pcs = arrayfun(@(i) {'Electrolyte', 'pcs', i}, (1 : model.Electrolyte.nlogsp)', 'uniformoutput', false);

            elyte = 'Electrolyte';
            indsolid = model.(elyte).indsolidsp(1);
            
            p = {qp{:}                             , ...
                 pcs{:}                            , ...
                 {'Anode', 'E'}                    , ...
                 {'Electrolyte', 'phi'}            , ...
                 {'Anode'      , 'volumeFraction'} , ...
                 {'Cathode'    , 'volumeFraction'} , ...
                 {'Electrolyte', 'volumeFraction'} , ...
                 {'Cathode'    , 'E'}};

        end
        
        function [state, report] = stepFunction(model, state, state0, dt, drivingForces, linsolver, nonlinsolver, iteration, varargin)
            
            [state, report] = stepFunction@SeaWaterBattery(model, state, state0, dt, drivingForces, linsolver, nonlinsolver, iteration, varargin{:});
            
        end
        
        function cleanState = addStaticVariables(model, cleanState, state)
            
            cleanState = addStaticVariables@MagnesiumBattery(model, cleanState, state);
           
            elyte = 'Electrolyte';
            indsolid = model.(elyte).indsolidsp(1);
            
            cleanState.(elyte).cs{indsolid} = state.(elyte).cs{indsolid};
            cleanState.(elyte).nucleation   = state.(elyte).nucleation;
            
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
