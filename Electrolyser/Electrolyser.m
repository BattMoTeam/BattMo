classdef Electrolyser < BaseModel
    
    properties
        
        con = PhysicalConstants();

        % Components
        IonomerMembrane
        HydrogenEvolutionElectrode
        OxygenEvolutionElectrode
        
        couplingTerms
        couplingNames

        primaryVarNames
        funcCallList
        
    end
    
    methods
        
        function model = Electrolyser(paramobj)

            model = model@BaseModel();

            model.G = paramobj.G;
            
            model.HydrogenEvolutionElectrode = EvolutionElectrode(paramobj.HydrogenEvolutionElectrode);
            model.OxygenEvolutionElectrode   = EvolutionElectrode(paramobj.OxygenEvolutionElectrode);
            model.IonomerMembrane            = IonomerMembrane(paramobj.IonomerMembrane);
            
        end
        
        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@BaseModel(model);

            varnames = {'T'};

            model = model.registerVarNames(varnames);

            inm = 'IonomerMembrane';
            her = 'HydrogenEvolutionElectrode';
            oer = 'OxygenEvolutionElectrode';

            ctl = 'CatalystLayer';
            exl = 'ExchangeLayer';
            ptl = 'PorousTransportLayer';
            
            fn = @Electrolyser.dispatchTemperature;
            inputvarnames = {'T'};
            model = model.registerPropFunction({{her, 'T'}, fn, inputvarnames});            
            model = model.registerPropFunction({{oer, 'T'}, fn, inputvarnames});
            model = model.registerPropFunction({{inm, 'T'}, fn, inputvarnames});
            
            fn = @Electrolyser.updateC;

            eldes = {her, oer};
            layers = {ctl, exl};
            
            fn = @Electrolyser.updateIonomerSources;
            inputvarnames = {};
            for ielde = 1 : numel(eldes)
                elde = eldes{ielde};
                inputvarnames{end + 1} = {elde, ctl, 'elyteReactionRate'};
                inputvarnames{end + 1} = {elde, ctl, 'inmrReactionRate'};
                inputvarnames{end + 1} = {elde, exl, 'OHexchangeRate'};
                inputvarnames{end + 1} = {elde, exl, 'H2OexchangeRate'};
            end
            model = model.registerPropFunction({{inm, 'H2OSource'}, fn, inputvarnames});
            model = model.registerPropFunction({{inm, 'OHSource'}, fn, inputvarnames});


            fn = @Electrolyser.dispatchIonomerToReactionLayers;
            inputvarnames = {{inm, 'phi'}, ...
                             {inm, 'cOH'}, ...
                             {inm, 'H2Oa'}};
            for ielde = 1 : numel(eldes)
                for ilayer = 1 : numel(layers)
                    elde = eldes{ielde};
                    layer = layers{ilayer};
                    model = model.registerPropFunction({{elde, layer, 'cOHinmr'}, fn, inputvarnames});
                    model = model.registerPropFunction({{elde, layer, 'phiInmr'}, fn, inputvarnames});
                    model = model.registerPropFunction({{elde, layer, 'H2OaInmr'}, fn, inputvarnames});
                end
            end
            
            model = model.registerStaticVarNames({{inm, 'jBcSource'}, ...
                                                  'T'});
            
        end


        function model = validateModel(model, varargin)

            model = validateModel@BaseModel(model, varargin{:});
            cgt = ComputationalGraphTool(model);

            model.primaryVarNames = cgt.getPrimaryVariables();
            model.funcCallList = cgt.setOrderedFunctionCallList();
            
        end
        
        function initstate = setupInitialState(model)
            
        end

        function state = dispatchTemperature(model, state)
            
            oer = 'OxygenEvolutionElectrode';
            her = 'HydrogenEvolutionElectrode';
            inm = 'IonomerMembrane';
            
            state.(oer).T = state.T(model.(oer).G.cellmap);
            state.(her).T = state.T(model.(her).G.cellmap);
            state.(inm).T = state.T(model.(inm).G.cellmap);
            
        end
        
        function state = dispatchIonomerToReactionLayers(model, state)

            inm = 'IonomerMembrane';
            her = 'HydrogenEvolutionElectrode';
            oer = 'OxygenEvolutionElectrode';

            ctl = 'CatalystLayer';
            exl = 'ExchangeLayer';

            eldes = {her, oer};

            for ielde = 1 : numel(eldes)
                
                elde = eldes{ielde};
                coupterms = model.couplingTerms;
                coupnames = model.couplingNames;
                coupterm  = getCoupTerm(coupterms, {inm, elde}, coupnames);
                coupcells = coupterm.couplingcells;

                for ilayer = 1 : numel(layers)
                    layer = layers{ilayer};
                    state.(elde).(layer).H2OaInmr(coupcells(:, 2)) = state.(inm).H2Oa(coupcells(:, 1));
                    state.(elde).(layer).cOHinmr(coupcells(:, 2))  = state.(inm).cOH(coupcells(:, 1));
                    state.(elde).(layer).phiInmr(coupcells(:, 2))  = state.(inm).phi(coupcells(:, 1));
                end
                
            end
            
        end

        function state = updateIonomerSources(model, state)
            
            inm = 'IonomerMembrane';
            her = 'HydrogenEvolutionElectrode';
            oer = 'OxygenEvolutionElectrode';

            ctl = 'CatalystLayer';
            exl = 'ExchangeLayer';

            % initialize sources (done in a way that takes care of AD, meaning that OHSource and H2OSource are AD
            % variables is state.(inm).phi is)
            OHSource = 0*state.(inm).phi;
            H2OSource = 0*state.(inm).phi;

            eldes = {her, oer};

            for ielde = 1 : numel(eldes)
                
                elde = eldes{ielde};
                coupterms = model.couplingTerms;
                coupnames = model.couplingNames;
                coupterm  = getCoupTerm(coupterms, {inm, elde}, coupnames);
                coupcells = coupterm.couplingcells;

                reacR    = state.(elde).(ctl).inmrReactionRate(coupcells(:, 2));
                OHexchR  = state.(elde).(exl).OHexchangeRate(coupcells(:, 2));
                H2OexchR = state.(elde).(exl).H2OexchangeRate(coupcells(:, 2));
                
                leps  = model.(inm).liquidVolumeFraction(coupcells(:, 2));
                
                OHSource(coupcells(:, 1))  = -2.*reacR./(n.*F) - OHexchR.*leps;
                H2OSource(coupcells(:, 1)) = -reacR./(n.*F) - H2OexchR;
                
            end
            
            state.(inm).OHSource  = OHSource;
            state.(inm).H2OSource = H2OSource;
            
        end

        function primaryvarnames = getPrimaryVariables(model)
            primaryvarnames = model.primaryVarNames;
        end


        function cleanState = addStaticVariables(model, cleanState, state)

            cleanState = addStaticVariables@BaseModel(model, cleanState, state);
            cleanState.T = state.T;            
            
        end

        
        function forces = getValidDrivingForces(model)
            
            forces = getValidDrivingForces@PhysicalModel(model);
            
        end

        function [problem, state] = getEquations(model, state0, state,dt, drivingForces, varargin)
            
            opts = struct('ResOnly', false, 'iteration', 0); 
            opts = merge_options(opts, varargin{:});
            
            time = state0.time + dt;
            if(not(opts.ResOnly))
                state = model.initStateAD(state);
            end

            funcCallList = model.funcCallList;

            for ifunc = 1 : numel(funcCallList)
                eval(funcCallList{ifunc});
            end
            
            %% Set up the governing equations
            
            eqs = {};
            eqs{end + 1} = state.IonomerMembrane.chargeCons;
            eqs{end + 1} = state.IonomerMembrane.H2OmassCons;
            eqs{end + 1} = state.HydrogenEvolutionElectrode.PorousTransportLayer.Boundary.gasStateEquation;
            eqs{end + 1} = state.HydrogenEvolutionElectrode.PorousTransportLayer.Boundary.liquidStateEquation;
            eqs{end + 1} = state.HydrogenEvolutionElectrode.PorousTransportLayer.chargeCons;
            eqs{end + 1} = state.HydrogenEvolutionElectrode.PorousTransportLayer.compGasMassCons{1};
            eqs{end + 1} = state.HydrogenEvolutionElectrode.PorousTransportLayer.compGasMassCons{2};
            eqs{end + 1} = state.HydrogenEvolutionElectrode.PorousTransportLayer.liquidMassCons;
            eqs{end + 1} = state.HydrogenEvolutionElectrode.PorousTransportLayer.OHMassCons;
            eqs{end + 1} = state.HydrogenEvolutionElectrode.PorousTransportLayer.liquidStateEquation;
            eqs{end + 1} = state.OxygenEvolutionElectrode.PorousTransportLayer.Boundary.gasStateEquation;
            eqs{end + 1} = state.OxygenEvolutionElectrode.PorousTransportLayer.Boundary.liquidStateEquation;
            eqs{end + 1} = state.OxygenEvolutionElectrode.PorousTransportLayer.chargeCons;
            eqs{end + 1} = state.OxygenEvolutionElectrode.PorousTransportLayer.compGasMassCons{1};
            eqs{end + 1} = state.OxygenEvolutionElectrode.PorousTransportLayer.compGasMassCons{2};
            eqs{end + 1} = state.OxygenEvolutionElectrode.PorousTransportLayer.liquidMassCons;
            eqs{end + 1} = state.OxygenEvolutionElectrode.PorousTransportLayer.OHMassCons;
            eqs{end + 1} = state.OxygenEvolutionElectrode.PorousTransportLayer.liquidStateEquation;

            neq = numel(eqs);
            
            types = repmat({'cell'}, 1, neq);
            names = arrayfun(@(ind) sprintf('eqs_%d', ind), (1 : neq)', 'uniformoutput', false);
            
            primaryVars = model.getPrimaryVariables();

            %% setup LinearizedProblem that can be processed by MRST Newton API
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
        
        end
        
    end
    
end


