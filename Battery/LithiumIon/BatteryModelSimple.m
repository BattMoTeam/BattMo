classdef BatteryModelSimple < PhysicalModel

    properties

        con = physicalConstants();

        % Coupling structures
        couplingnames
        couplingTerms

        % Temperature and SOC
        % for the moment here, for convenience. Will be moved
        T
        SOC

        % Input current
        J
        % Voltage cut
        Ucut

        % Submodels
        elyte
        ne
        pe
        ccpe
        ccne
        sep
        
    end

    methods

        function model = BatteryModelSimple(params)

            model = model@PhysicalModel([]);
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks',false);

            %% Setup the model using the input parameters
            model.T     = params.T;
            model.SOC   = params.SOC;
            model.J     = params.J;
            model.Ucut  = params.Ucut;

            model.G = params.G;

            model.elyte = params.elyte;
            model.ne    = params.ne;
            model.pe    = params.pe;
            model.ccpe  = params.ccpe;
            model.ccne  = params.ccne;
            model.sep   = params.sep;
            
            model = model.setElytePorosity();
            
            %% Setup the couplings using the input parameters
            coupTerms = {};
            coupTerms{end + 1} = params.coupTermNeElyte;
            coupTerms{end + 1} = params.coupTermPeElyte;
            coupTerms{end + 1} = params.coupTermCcneNe;
            coupTerms{end + 1} = params.coupTermCcpePe;
            coupTerms{end + 1} = params.coupTermCcneBc;
            coupTerms{end + 1} = params.coupTermCcpeBc;
            model.couplingTerms = coupTerms;
            model.couplingnames = cellfun(@(x) x.name, coupTerms, 'uniformoutput', false);
            
        end
        
        function [problem, state] = getEquations(model, state0, state,dt, drivingForces, varargin)
            
            time = state0.time+dt;
            state=model.initStateAD(state);
            
            %% for now temperature and SOC are kept constant
            nc = model.G.cells.num;
            state.T   = model.T*ones(nc, 1);
            state.SOC = model.SOC*ones(nc, 1);
            
            state = model.dispatchValues(state);
            state = model.updatePhiElyte(state);
            
            %% Update Source and BC term variables
            
            names={{'pe','am'},{'pe','am'}};
            for i=1:numel(names)
                submodel = model.getSubmodel(names{i});
                val = submodel.updateQuantities(model.getProp(state,names{i}));
                state = model.setProp(state,names{i},val);
            end
            
            state = setupBCSources(model, state);
            
            %% Update Reaction Coupling variables
            
            names={{'ne'},{'pe'}};
            for i=1:numel(names)
                submodel=model.getSubmodel(names{i});
                val = submodel.updateReactionRate(model.getProp(state,names{i}));
                state = model.setProp(state,names{i},val);
            end
            
            state = setupExchanges(model, state);
            
            %% Update Fluxes variables
            
            names={{'elyte'},{'ne'},{'pe'}};
            for i=1:numel(names)
                submodel=model.getSubmodel(names{i});
                val = submodel.updateQuantities(model.getProp(state,names{i}));
                state = model.setProp(state,names{i},val);
            end
            
            names={{'ccpe'},{'ccne'}};
            for i=1:numel(names)
                submodel=model.getSubmodel(names{i});
                val = submodel.updateChargeCont(model.getProp(state,names{i}));
                state = model.setProp(state,names{i},val);
            end
            
            %% Set up the governing equations
            
            %% Mass and charge conservation for the electorlyte and the electrode
            
            % Accumulation terms for the mass conservation equtions
            cdotLi=struct();
            cdotLi.elyte = (state.elyte.cs{1} - state0.elyte.cs{1})/dt;
            cdotLi.ne    = (state.ne.am.Li - state0.ne.am.Li)/dt;
            cdotLi.pe    = (state.pe.am.Li - state0.pe.am.Li)/dt;
            
            names={'elyte','ne','pe'};
            eqs={};
            
            for i = 1 : numel(names)
                
                submodel = model.getSubmodel({names{i}});

                %% probably only be done on the submodel
                source = model.getProp(state,{names{i},'LiSource'});
                flux = model.getProp(state,{names{i},'LiFlux'});

                %% could use submodel
                div = submodel.operators.Div(flux)./submodel.G.cells.volumes;

                %% HAC
                if(strcmp(names{i},'elyte'))
                    cepsdot = submodel.eps.*cdotLi.(names{i});
                else
                    cepsdot = submodel.am.eps.*cdotLi.(names{i});
                end
                %% Li conservation
                eqs{end+1} = -div + source - cepsdot;
                
                %¤ charge continutity
                %% should probably be done on the sub model
                eqs{end+1} = model.getProp(state,{names{i},'chargeCont'});
            end
            
            %% charge conservation for the current collectors
            
            names = {'ccne','ccpe'};
            for i = 1 : numel(names)
                eqs{end+1} = model.getProps(state, {names{i}, 'chargeCont'});
            end
            
            %% setup control equation (fixed total current at ccpe)
            
            src = drivingForces.src(time);
            coupterm = model.getCoupTerm('bc-ccpe');
            faces = coupterm.couplingfaces;
            bcval = state.ccpe.E;
            ccpe_sigmaeff = model.ccpe.sigmaeff;
            [tccpe, cells] = model.ccpe.operators.harmFaceBC(ccpe_sigmaeff, faces);
            control = src - sum(tccpe.*(bcval - state.ccpe.phi(cells)));
            
            eqs{end+1} = control;
            
            
            %% Give type and names to equations and names of the primary variables (for book-keeping)
            
            types={'cell','cell','cell','cell',...
                   'cell','cell','cell','cell','cell'};
            names = {'elyte_Li_massCont', ...
                     'elyte_chargeCont' , ...
                     'ne_Li_massCont'   , ...
                     'ne_e_chargeCont'  , ...
                     'pe_Li_massCont'   , ...
                     'pe_e_chargeCont'  , ...
                     'ccne_e_chargeCont', ...
                     'ccpe_e_chargeCont', ...
                     'control'};
            primaryVars = model.getPrimaryVariables();

            %% setup LinearizedProblem that can be processed by MRST Newton API
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
            
        end
        
        function state = dispatchValues(model, state)

            T = state.T;
            SOC = state.SOC;

            elyte = model.elyte;
            G = elyte.G;
            Telyte = T(G.mappings.cellmap);

            ne = model.ne;
            G = ne.G;
            Tne = T(G.mappings.cellmap);
            SOCne = SOC(G.mappings.cellmap);
            
            pe = model.pe;
            G = pe.G;
            Tpe = T(G.mappings.cellmap);
            SOCpe = SOC(G.mappings.cellmap);
            
            ccpe = model.ccpe;
            G = ccpe.G;
            Tccpe = T(G.mappings.cellmap);

            ccne = model.ccne;
            G = ccne.G;
            Tccne = T(G.mappings.cellmap);

            state.elyte.T = Telyte;
            state.ne.T    = Tne;
            state.ne.am.T = Tne;
            state.pe.T    = Tpe;
            state.pe.am.T = Tpe;
            state.ccpe.T  = Tccpe;
            state.ccne.T  = Tccne;
            
            state.ne.SOC    = SOCne;
            state.ne.am.SOC = SOCne;
            state.pe.SOC    = SOCpe;
            state.pe.am.SOC = SOCpe;
            
        end

        function state = updatePhiElyte(model, state)

            phielyte = state.elyte.phi;
            elyte = model.elyte;

            elytecells = zeros(model.G.cells.num, 1);
            elytecells(elyte.G.mappings.cellmap) = (1 : elyte.G.cells.num)';

            ne = model.ne;
            phielyte_ne = phielyte(elytecells(ne.G.mappings.cellmap));

            pe = model.pe;
            phielyte_pe = phielyte(elytecells(pe.G.mappings.cellmap));

            state.ne.phielyte = phielyte_ne;
            state.pe.phielyte = phielyte_pe;

        end

        
        function model = setElytePorosity(model)

            elyte = model.getSubmodel({'elyte'});
            ne = model.getSubmodel({'ne'});
            pe = model.getSubmodel({'pe'});
            sep = model.getSubmodel({'sep'});

            elytecells = zeros(model.G.cells.num, 1);
            elytecells(elyte.G.mappings.cellmap) = (1 : elyte.G.cells.num)';

            elyte.eps = NaN(elyte.G.cells.num, 1);
            elyte.eps(elytecells(ne.G.mappings.cellmap)) = ne.void;
            elyte.eps(elytecells(pe.G.mappings.cellmap)) = pe.void;
            elyte.eps(elytecells(sep.G.mappings.cellmap)) = sep.void;

            model.elyte =elyte;

        end
        
        function initstate = setupInitialState(model)
            
            % Setup initial state

            nc = model.G.cells.num;

            SOC = model.SOC;
            T   = model.T;
            
            initstate.T   =  T*ones(nc, 1);
            initstate.SOC =  SOC*ones(nc, 1);

            initstate = model.dispatchValues(initstate);
            
            elyte = model.elyte;
            ne    = model.ne;
            pe    = model.pe;
            ccne  = model.ccne;
            ccpe  = model.ccpe;

            ne_am = ne.am;
            pe_am = pe.am;

            %% setup initial ne state

            m = (1 ./ (ne_am.theta100 - ne_am.theta0));
            b = -m .* ne_am.theta0;
            theta = (SOC - b) ./ m;
            c = theta .* ne_am.Li.cmax;
            c = c*ones(ne.G.cells.num, 1);

            initstate.ne.am.Li = c;
            initstate.ne.am = ne_am.updateQuantities(initstate.ne.am);

            OCP = initstate.ne.am.OCP;
            initstate.ne.am.phi = OCP;

            %% setup initial pe state

            m = (1 ./ (pe_am.theta100 - pe_am.theta0));
            b = -m .* pe_am.theta0;
            theta = (SOC - b) ./ m;
            c = theta .* pe_am.Li.cmax;
            c = c*ones(ne.G.cells.num, 1);

            initstate.pe.am.Li = c;
            initstate.pe.am = pe_am.updateQuantities(initstate.pe.am);

            OCP = initstate.pe.am.OCP;
            initstate.pe.am.phi = OCP;

            %% setup initial elyte state

            initstate.elyte.phi = zeros(elyte.G.cells.num, 1);
            cs=cell(2,1);
            initstate.elyte.cs=cs;
            initstate.elyte.cs{1} = 1000*ones(elyte.G.cells.num, 1);

            %% setup initial Current collectors state
            
            OCP = initstate.ne.am.OCP;
            OCP = OCP(1) .* ones(ccne.G.cells.num, 1);
            initstate.ccne.phi = OCP;

            OCP = initstate.pe.am.OCP;
            OCP = OCP(1) .* ones(ccpe.G.cells.num, 1);
            initstate.ccpe.phi = OCP;
            
            initstate.ccpe.E = OCP(1);
            
        end
        
        function coupterm = getCoupTerm(model, coupname)
            coupnames = model.couplingnames;
            
            [isok, ind] = ismember(coupname, coupnames);
            assert(isok, 'name of coupling term is not recognized.');
            
            coupterm = model.couplingTerms{ind};
            
        end
        
        function state = initStateAD(model,state)
            adbackend = model.AutoDiffBackend();
            [state.elyte.cs{1},...
             state.elyte.phi,...   
             state.ne.am.Li,...    
             state.ne.am.phi,...   
             state.pe.am.Li,...    
             state.pe.am.phi,...   
             state.ccne.phi,...    
             state.ccpe.phi,...    
             state.ccpe.E]=....
                adbackend.initVariablesAD(...
                    state.elyte.cs{1},...
                    state.elyte.phi,...   
                    state.ne.am.Li,...    
                    state.ne.am.phi,...   
                    state.pe.am.Li,...    
                    state.pe.am.phi,...   
                    state.ccne.phi,...    
                    state.ccpe.phi,...    
                    state.ccpe.E);       
        end
        
        function p = getPrimaryVariables(model)
            p ={{'elyte','cs'},...
                {'elyte','phi'},...   
                {'ne','am','Li'},...    
                {'ne','am','phi'},...   
                {'pe','am','Li'},...    
                {'pe','am','phi'},...   
                {'ccne','phi'},...    
                {'ccpe','phi'},...    
                {'ccpe','E'}
               };
        end
        
        function state = setProp(model,state,names,val)
            nname=numel(names);
            if(nname==1)
                state.(names{1})=val;
            elseif(nname==2)
                state.(names{1}).(names{2})=val;
            elseif(nname==3)
                state.(names{1}).(names{2}).(names{3}) = val;
            else
                error('not implmented')
            end
        end
        
        function var = getProp(model, state,names)
            var = state.(names{1});
            for i=2:numel(names)
                var = var.(names{i});
            end
        end
        
        function submod = getSubmodel(model, names)
            submod = model.(names{1});
            for i=2:numel(names)
                submod = submod.(names{i});
            end
        end

        function validforces = getValidDrivingForces(model)
            validforces=struct('src', [], 'stopFunction', []); 
        end
         
        
        function [state, report] = updateState(model,state, problem, dx, drivingForces)
            p = model.getPrimaryVariables();
            for i=2:numel(dx)
                val = model.getProps(state,p{i});
                val = val + dx{i};
                state = model.setProp(state,p{i},val);
            end
            %% not sure how to handle cells
            state.elyte.cs{1} =  state.elyte.cs{1} + dx{1};
            report = [];
        end
        
        
        function state = reduceState(model, state, removeContainers)
        % Reduce state to double (do not use property containers)
            state = value(state);
        end

    end

end
