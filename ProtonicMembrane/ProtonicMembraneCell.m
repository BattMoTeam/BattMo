classdef ProtonicMembraneCell < BaseModel
    
    properties
        
        % Temperature
        T
        % Structure with physical constants
        constants

        Anode
        Cathode
        Electrolyte
        
        couplingTerms
        couplingnames

        iBV
        
    end
    
    methods
        
        function model = ProtonicMembraneCell(paramobj)

            model = model@BaseModel();

            fdnames = {'couplingTerms'};
            model = dispatchParams(model, paramobj, fdnames);

            model.Anode       = ProtonicMembraneElectrode(paramobj.Anode);
            model.Cathode     = ProtonicMembraneElectrode(paramobj.Cathode);
            model.Electrolyte = ProtonicMembraneElectrolyte(paramobj.Electrolyte);
            
            % setup couplingNames
            model.couplingnames = cellfun(@(x) x.name, model.couplingTerms, 'uniformoutput', false);

            model = model.setupOperators();
            
            % setup standard physical constants
            model.constants = PhysicalConstants();

        end
        
        function model = registerVarAndPropfuncNames(model)
            
            an    = 'Anode';
            ct    = 'Cathode';
            elyte = 'Electrolyte';
            
            model = registerVarAndPropfuncNames@BaseModel(model);
            
            varnames = {};
            % control equations
            varnames{end + 1} = 'controlEqs';
            % control equations
            varnames{end + 1} = 'boundaryEqs';
            
            model = model.registerVarNames(varnames);
        
            fn = @ProtonicMembraneCell.dispatchE;
            inputnames = {{elyte, 'E'}};
            model = model.registerPropFunction({{an, 'E'}, fn, inputnames});
            model = model.registerPropFunction({{ct, 'E'}, fn, inputnames});
            
            fn = @ProtonicMembraneCell.setupHpSources;
            inputnames = {{an, 'phiElectrolyte'}, {ct, 'phiElectrolyte'}, {'elyte', 'phi'}};
            model = model.registerPropFunction({{elyte, 'sourceHp'}, fn, inputnames});
            
            fn = @ProtonicMembraneCell.setupElSources;
            inputnames = {{an, 'E'}, {an, 'phiElectrolyte'}, ...
                          {ct, 'E'}, {ct, 'phiElectrolyte'}, ...
                          {elyte, 'E'}, {elyte, 'phi'}, ...
                          {elyte, 'sigmaEl'}};
            model = model.registerPropFunction({{elyte, 'sourceEl'}, fn, inputnames});

            fn = @ProtonicMembraneCell.updateControlEqs;
            inputnames = {{an, 'iBV'}, {ct, 'E'}};
            model = model.registerPropFunction({'controlEqs', fn, inputnames});
            
            fn = @ProtonicMembraneCell.updateBoundaryEqs;
            inputnames = {{an, 'iBV'}, {ct, 'iBV'}, {elyte, 'sourceHp'}};
            model = model.registerPropFunction({'boundaryEqs', fn, inputnames});

            
        end

        function model = setupOperators(model)

        % this function is  hacky at the moment!! We will make it better later

            an    = 'Anode';
            ct    = 'Cathode';
            elyte = 'Electrolyte';
            
            G = model.(elyte).G;
            
            couplingterms = model.couplingTerms;
            coupnames     = model.couplingnames;

            % Anode part
            coupterm = getCoupTerm(couplingterms, 'Anode-Electrolyte', coupnames);
            d = G.cells.centroids(coupterm.couplingcells(:, 2), :) -  G.faces.centroids(coupterm.couplingfaces(:, 2), :);
            d = sqrt(sum(d.^2, 2));
            anode_T = 1./d; 
            
            % Cathode part
            coupterm = getCoupTerm(couplingterms, 'Cathode-Electrolyte', coupnames);
            d = G.cells.centroids(coupterm.couplingcells(:, 2), :) - G.faces.centroids(coupterm.couplingfaces(:, 2), :);
            d = sqrt(sum(d.^2, 2));
            cathode_T = 1./d; 
            
            model.operators.anode_T   = anode_T;
            model.operators.cathode_T = cathode_T;
            
        end
        
        function state = setupHpSources(model, state)

            an    = 'Anode';
            ct    = 'Cathode';
            elyte = 'Electrolyte';
            
            coupterms = model.couplingTerms;
            coupnames = model.couplingNames;
            sigmaHp   = model.sigmaHp;
            
            an_phi   = state.(an).phiElectrolyte;
            ct_phi   = state.(ct).phiElectrolyte;
            phi      = state.(elyte).phi;
            
            sourceHp = 0*phi; % update sourceHp as AD
            
            % Anode part
            coupterm = getCoupTerm(coupterms, 'Anode-Electrolyte', coupnames);
            cc = coupterm.couplingcells(:, 2);
            sourceHp(cc) = op.anode_T(cc).*sigmaHp.*(an_phi(coupterm.couplingcells(:, 1)) - phi(cc));
            
            % Cathode part
            coupterm = getCoupTerm(coupterms, 'Cathode-Electrolyte', coupnames);
            cc = coupterm.couplingcells(:, 2);
            sourceHp(cc) = op.cathode_T(cc).*sigmaHp.*(ct_phi(coupterm.couplingcells(:, 1)) - phi(cc));
            
        end



        function state = setupElSources(model, state)
            
            an    = 'Anode';
            ct    = 'Cathode';
            elyte = 'Electrolyte';
            
            op = model.operators
            
            an_E    = state.(an).E;
            ct_E    = state.(ct).E;
            E       = state.(elyte).E;
            an_phi  = state.(an).phiElectrolyte;
            ct_phi  = state.(ct).phiElectrolyte;
            phi     = state.(elyte).phi;
            sigmaEl = state.(elyte).sigmaEl;
            
            sourceEl = 0*state.(elyte).phi; % initialize AD for sourceEl
            
            % Anode part
            coupterm = getCoupTerm(couplingterms, 'Anode-Electrolyte', coupnames);
            eldecc = coupterm.couplingcells(:, 1);
            cc = coupterm.couplingcells(:, 2);
            sourceEl(cc) = op.anode_T(cc).*sigmaEl(cc).*(an_phi(eldecc) + an_E(eldecc) - (E(cc) + phi(cc)));
            
            % Cathode part
            coupterm = getCoupTerm(couplingterms, 'cathode-Electrolyte', coupnames);
            eldecc = coupterm.couplingcells(:, 1);
            cc = coupterm.couplingcells(:, 2);
            sourceEl(cc) = op.anode_T(cc).*sigmaEl(cc).*(an_phi(eldecc) + an_E(eldecc) - (E(cc) + phi(cc)));

        end

        
        function state = updateControlEqs(model, state)

            an    = 'Anode';
            ct    = 'Cathode';
            elyte = 'Electrolyte';
            
            eqs1 = state.(an).iBv - model.iBV;
            eqs2 = state.(ct).E;
            
            state.controlEqs = [eqs1; eqs2];
            
        end

        
        function state = updateBoundaryEqs(model, state)

            an    = 'Anode';
            ct    = 'Cathode';
            elyte = 'Electrolyte';
            
            coupterms = model.couplingTerms;
            coupnames = model.couplingNames;
            sigmaHp   = model.sigmaHp;
            
            an_phi   = state.(an).phiElectrolyte;
            phi      = state.(elyte).phi;
            sourceHp = state.(elyte).sourceHp;
            
            an    = 'Anode';
            ct    = 'Cathode';
            elyte = 'Electrolyte';
            
            an_iBV   = state.(an).iBV;
            ct_iBV   = state.(ct).iBV;
            sourceHp = state.(elyte).sourceHp;
            
            couplingterms = model.couplingTerms;
            coupnames     = model.couplingNames;
            
            % Anode part
            coupterm = getCoupTerm(couplingterms, 'Anode-Electrolyte', coupnames);
            eqs1 = sourceHp(coupterm.couplingcells(:, 2)) - an_iBV(coupterm.couplingcells(:, 1));
            
            % Cathode part
            coupterm = getCoupTerm(couplingterms, 'Cathode-Electrolyte', coupnames);
            eqs2 = sourceHp(coupterm.couplingcells(:, 2)) - ct_iBV(coupterm.couplingcells(:, 1));

            state.boundaryEqs = [eqs1; eqs2];
        
        end

        function primaryvarnames = getPrimaryVariables(model)
            
            an    = 'Anode';
            ct    = 'Cathode';
            elyte = 'Electrolyte';
            
            primaryvarnames = {{elyte, 'phi'} , ...
                               {elyte, 'E'}   , ...
                               {an, 'E'}      , ...
                               {an, 'phiElectrolyte'}    , ...
                               {ct, 'E'}
                              };
            
        end

        
        function cleanState = addStaticVariables(model, cleanState, state)

            cleanState = addStaticVariables@BaseModel(model, cleanState, state);
            
            cleanState.Cathode.phiElectrolyte = 0;
            cleanState.Cathode.Eocp = 0;
            cleanState.Anode.Eocp = 0;
            
        end
        
        function [problem, state] = getEquations(model, state0, state,dt, drivingForces, varargin)

            opts = struct('ResOnly', false, 'iteration', 0,'reverseMode',false); 
            opts = merge_options(opts, varargin{:});
            
            time = state0.time + dt;
            if(not(opts.ResOnly) && not(opts.reverseMode))
                state = model.initStateAD(state);
            elseif(opts.reverseMode)
                disp('No AD initatlization in equation old style')
                state0 = model.initStateAD(state0);
            else
                assert(opts.ResOnly);
            end

            state.(elyte) = model.elyte.updateElConductivity(state.(elyte));
            state.(elyte) = model.elyte.updateElFlux(state.(elyte));
            state.(elyte) = model.elyte.updateHpFlux(state.(elyte));
            state         = model.dispatchE(state);
            state.(an)    = model.an.updateButtlerVolmerRate(state.(an));
            state         = model.setupHpSources(state);
            state.(elyte) = model.elyte.updateChargeConsHp(state.(elyte));
            state         = model.dispatchE(state);
            state         = model.setupElSources(state);
            state.(elyte) = model.elyte.updateChargeConsEl(state.(elyte));
            state.(ct)    = model.ct.updateButtlerVolmerRate(state.(ct));
            state         = model.updateBoundaryEqs(state);
            state         = model.updateControlEqs(state);
            
            eqs = {};
            eqs{end + 1} = state.(elyte).chargeConsHp;
            eqs{end + 1} = state.(elyte).chargeConsEl;
            eqs{end + 1} = state.boundaryEqs;
            eqs{end + 1} = state.controlEqs;
            
            names = {'Hp_charge_cons', 'El_charge_cons', 'boundary_eqs', 'control_eqs'};
            types = repmat({'cells'}, 1, numel(names));
            
            primaryVars = model.getPrimaryVariables();
            
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
            
        end
        
        
    end
    
end
