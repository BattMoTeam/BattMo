classdef MagnesiumBattery < SeaWaterBattery

    properties
        
        doposqpcapping = true % Cap quasi-particle concentrations in the Newton update (only for the "always positive" quasiparticle)
        printVoltage = true   % print voltage reached at each converged step

    end
    
    
    methods

        function model = MagnesiumBattery(inputparams)
            
            model = model@SeaWaterBattery(inputparams);

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
            
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
                    
            opts = struct('ResOnly', false, 'iteration', 0);
            opts = merge_options(opts, varargin{:});

            elyte = 'Electrolyte';
            anam  = 'AnodeActiveMaterial';
            ctam  = 'CathodeActiveMaterial';
            ct    = 'Cathode';
            an    = 'Anode';

            time = state0.time + dt;
            if(not(opts.ResOnly))
                state = model.initStateAD(state);
            end

            funcCallList = model.funcCallList;

            for ifunc = 1 : numel(funcCallList)
                eval(funcCallList{ifunc});
            end
                        
            anvols    = model.(an).G.cells.volumes;
            ctvols    = model.(ct).G.cells.volumes;
            elytevols = model.(elyte).G.cells.volumes;
            F         = model.con.F;
            nqp       = model.(elyte).nqp;
            is_prep   = model.include_precipitation;
            
            eqs = {};
            for ind = 1 : nqp
                eqs{end + 1} = 1e-2./elytevols.*state.(elyte).qpMassCons{ind};
            end
            for ind = 1 : nqp
                eqs{end + 1} = 1e-2*state.(elyte).atomicMassCons{ind};
            end
            eqs{end + 1} = 1/F*state.(elyte).chargeCons;
            eqs{end + 1} = state.(ct).galvanostatic;
            eqs{end + 1} = state.(an).galvanostatic;

            if is_prep
                % Add the extra equations in case where precipitation is included.
                eqs{end + 1} = 1./anvols.*state.(an).massCons;
                eqs{end + 1} = 1./ctvols.*state.(ct).massCons;
                eqs{end + 1} = 1e-2./elytevols.*state.(elyte).dischargeMassCons;
                eqs{end + 1} = 1./elytevols.*state.(elyte).nucleationEquation;
                eqs{end + 1} = state.(elyte).volumeFractionEquation;
            end
            
            %% setup equation names (just for printed output)
            names = {};
            for ind = 1 : nqp
                names{end + 1} = sprintf('qp_%d_masscons', ind);
            end
            for ind = 1 : nqp
                names{end + 1} = sprintf('qp_%d_atommasscons', ind);
            end
            names = {names{:}, 'elyte_chargecons', 'ct_galvanostatic', 'an_galvanostatic'};

            if is_prep
                names = {names{:}, 'anode_masscons', 'cathode_masscons', 'discharge_masscons', 'nucleation', 'electrolyte_vf'};
            end
            
            types = repmat({'cell'}, 1, numel(names));

            primaryVars = model.getPrimaryVariableNames();

            %% setup LinearizedProblem that can be processed by MRST Newton API
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

        end

        %% Assembly functions
        
        function state = updateQuasiParticleSource(model, state)
        % Update source term for all the quasi particles

            coupDict  = model.couplingCellDict;
            qpdict    = model.Electrolyte.qpdict;
            nqp       = model.Electrolyte.nqp;
            nc        = model.Electrolyte.G.cells.num;
            vols      = model.Electrolyte.G.cells.volumes;
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
            vols = model.(ct).G.cells.volumes;

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

