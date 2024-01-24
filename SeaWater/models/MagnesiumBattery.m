classdef MagnesiumBattery < SeaWaterBattery

    properties
        
        doposqpcapping = true % Cap quasi-particle concentrations in the Newton update (only for the "always positive" quasiparticle)
        printVoltage   = true   % print voltage reached at each converged step

    end
    
    
    methods

        function model = MagnesiumBattery(inputparams)
            
            model = model@SeaWaterBattery(inputparams);
            model = model.setupForSimulation();
            
        end
        
        function model = setupElectrolyte(model, inputparams)

            if inputparams.include_precipitation
                model.Electrolyte = MagnesiumElectrolyte(inputparams.Electrolyte);
            else
                model.Electrolyte = MagnesiumElectrolyteNoPrecipitation(inputparams.Electrolyte);
            end
                
        end

        function model = setupAnode(model, inputparams)
            
            model.Anode = MagnesiumElectrode(inputparams);
            
        end
        
        function model = setupAnodeActiveMaterial(model, inputparams)
            
            model.AnodeActiveMaterial = MagnesiumActiveMaterial(inputparams);
            
        end
            
        function model = registerVarAndPropfuncNames(model)
                    
            model = registerVarAndPropfuncNames@SeaWaterBattery(model);
            
            ct    = 'Cathode';
            ctam  = 'CathodeActiveMaterial';
            an    = 'Anode';
            anam  = 'AnodeActiveMaterial';
            elyte = 'Electrolyte';

            model = model.registerVarNames({{ct, 'galvanostatic'},
                                            {an, 'galvanostatic'}
                                           });

            % Update source term for all the quasi particles
            nqp = model.(elyte).nqp;
            fn = @() MagnesiumBattery.updateQuasiParticleSource;
            inputnames = {{ctam, 'R'} , ...
                          {anam, 'R'}};
            model = model.registerPropFunction({VarName({elyte}, 'qpSrcTerms', nqp), fn, inputnames});            

            fn = @() SeaWaterBattery.updateGalvanostatic;
            fn = {fn, @(propfunction) PropFunction.drivingForceFuncCallSetupFn(propfunction)};
            inputnames = {{ctam, 'R'}, {an, 'E'}};
            model = model.registerPropFunction({{ct, 'galvanostatic'}, fn, inputnames});
            model = model.registerPropFunction({{an, 'galvanostatic'}, fn, inputnames});
            
            eldes = {an, ct};
            for ielde = 1 : numel(eldes)
                elde = eldes{ielde};
                model = model.removeVarNames({{elde, 'chargeCons'}, ...
                                              {elde, 'eSource'}, ...
                                              {elde, 'conductivity'}, ...
                                              {elde, 'j'}, ...
                                              {elde, 'jBcSource'}});
            end

            model = model.removeVarName({ct, 'volumeFraction'});
            
        end

        function model = setupForSimulation(model)

            shortNames ={'Electrolyte'          , 'elyte';
                         'AnodeActiveMaterial'  , 'anam';
                         'CathodeActiveMaterial', 'ctam';
                         'Cathode'              , 'ct';
                         'Anode'                , 'an'};

            model = model.equipModelForComputation('shortNames', shortNames);

            elyte = 'Electrolyte';
            anam  = 'AnodeActiveMaterial';
            ctam  = 'CathodeActiveMaterial';
            ct    = 'Cathode';
            an    = 'Anode';

            anvols    = model.(an).G.getVolumes();
            ctvols    = model.(ct).G.getVolumes();
            elytevols = model.(elyte).G.getVolumes();
            F         = model.con.F;
            nqp       = model.(elyte).nqp;
            is_prep   = model.include_precipitation;

            scalings = {};
            for ind = 1 : nqp
                scalings = horzcat(scalings, ...
                                   {{{'elyte', 'qpMassCons', ind}, elytevols/1e-2}, ...
                                    {{'elyte', 'atomicMassCons', ind}, 1/1e-2}});
            end
            scalings = horzcat(scalings, ...
                               {{{'elyte', 'chargeCons'}, F}});

            if is_prep
                scalings = horzcat(scalings                                       , ...
                                   {{{an, 'massCons'}, anvols}                    , ...
                                    {{ct, 'massCons'}, ctvols}                    , ...
                                    {{elyte, 'dischargeMassCons'}, elytevols/1e-2}, ...
                                    {{elyte, 'nucleationEquation'}, elytevols}});
            end
            
        end

        
        function state = updateQuasiParticleSource(model, state)
        % Update source term for all the quasi particles

            coupDict  = model.couplingCellDict;
            qpdict    = model.Electrolyte.qpdict;
            nqp       = model.Electrolyte.nqp;
            nc        = model.Electrolyte.G.getNumberOfCells();
            vols      = model.Electrolyte.G.getVolumes();
            adbackend = model.AutoDiffBackend;
            
            cathR    = state.CathodeActiveMaterial.R;
            anR      = state.AnodeActiveMaterial.R;

            % initialize the source terms
            for ind = 1 : nqp
                qpSrcTerms{ind} = adbackend.convertToAD(zeros(nc, 1), cathR);
            end

            coupcells = coupDict('Cathode-Electrolyte');
            ind = qpdict('HOH');
            qpSrcTerms{ind}(coupcells(:, 2)) = 2*vols(coupcells(:, 2)).*cathR(coupcells(:, 1));
            
            coupcells = coupDict('Anode-Electrolyte');
            ind = qpdict('Mg');
            qpSrcTerms{ind}(coupcells(:, 2)) = vols(coupcells(:, 2)).*anR(coupcells(:, 1));
            
            state.Electrolyte.qpSrcTerms = qpSrcTerms;
            
        end

        function state = updateGalvanostatic(model, state, drivingForces)

            ctam = 'CathodeActiveMaterial';
            ct   = 'Cathode';
            an   = 'Anode';
            
            t = state.time;
            R = state.(ctam).R;

            I = drivingForces.src(t);
            F = model.con.F;
            vols = model.(ct).G.getVolumes();

            state.(ct).galvanostatic = I + sum(2*F*R.*vols);

            state.(an).galvanostatic = state.(an).E;
            
        end

        %% Specialized Newton API methods
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)

            [state, report] = updateAfterConvergence@SeaWaterBattery(model, state0, state, dt, drivingForces);
            
            if model.printVoltage
                fprintf('\n\n *** voltage : %g\n\n', state.Cathode.E);
            end
            
        end

        function [state, report] = updateState(model, state, problem, dx, drivingForces)

            [state, report] = updateState@SeaWaterBattery(model, state, problem, dx, drivingForces);

            doposqpcapping = model.doposqpcapping;
            verbose = model.newton_verbose;
            
            if doposqpcapping

                elyte = 'Electrolyte';
                
                ics = model.(elyte).qpdict('Mg');           
                if min(state.(elyte).qpepscs{ics}) < 0
                    if verbose
                        fprintf('\n** cap total Mg, min value %g\n\n', min(state.(elyte).qpepscs{ics}));
                    end
                    state.(elyte).qpepscs{ics} = max(0, state.(elyte).qpepscs{ics});
                end
                
                ics = model.(elyte).qpdict('Cl');           
                if min(state.(elyte).qpepscs{ics}) < 0
                    if verbose
                        fprintf('\n** cap total Cl, min value %g\n\n', min(state.(elyte).qpepscs{ics}));
                    end
                    state.(elyte).qpepscs{ics} = max(0, state.(elyte).qpepscs{ics});
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
