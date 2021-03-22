classdef BatteryModelSimple < PhysicalModel

    properties

        con = PhysicalConstants();

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
        Electrolyte
        NegativeElectrode
        PositiveElectrode
        PositiveCurrentCollector
        NegativeCurrentCollector
        sep
        
    end

    methods

        function model = BatteryModelSimple(params)

            model = model@PhysicalModel([]);
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks',true);

            %% Setup the model using the input parameters
            model.T     = params.T;
            model.SOC   = params.SOC;
            model.J     = params.J;
            model.Ucut  = params.Ucut;

            model.G = params.G;

            model.Electrolyte              = params.Electrolyte;
            model.NegativeElectrode        = params.NegativeElectrode;
            model.PositiveElectrode        = params.PositiveElectrode;
            model.PositiveCurrentCollector = params.PositiveCurrentCollector;
            model.NegativeCurrentCollector = params.NegativeCurrentCollector;
            model.sep                      = params.sep;
            
            model = model.setElectrolytePorosity();
            
            %% Setup the couplings using the input parameters
            coupTerms = {};
            coupTerms{end + 1} = params.coupTermNegativeElectrodeElectrolyte;
            coupTerms{end + 1} = params.coupTermPositiveElectrodeElectrolyte;
            coupTerms{end + 1} = params.coupTermNegativeCurrentCollectorNegativeElectrode;
            coupTerms{end + 1} = params.coupTermPositiveCurrentCollectorPositiveElectrode;
            coupTerms{end + 1} = params.coupTermNegativeCurrentCollectorBc;
            coupTerms{end + 1} = params.coupTermPositiveCurrentCollectorBc;
            model.couplingTerms = coupTerms;
            model.couplingnames = cellfun(@(x) x.name, coupTerms, 'uniformoutput', false);
            
        end
        
        function [problem, state] = getEquations(model, state0, state,dt, drivingForces, varargin)
            
            time = state0.time + dt;
            state = model.initStateAD(state);
            
            %% for now temperature and SOC are kept constant
            nc = model.G.cells.num;
            state.T   = model.T*ones(nc, 1);
            state.SOC = model.SOC*ones(nc, 1);
            
            state = model.dispatchValues(state);
            state = model.updatePhiElectrolyte(state);
            
            %% Update Source and BC term variables
            
            names={{'PositiveElectrode','ActiveMaterial'},{'PositiveElectrode','ActiveMaterial'}};
            for i=1:numel(names)
                submodel = model.getSubmodel(names{i});
                val = submodel.updateQuantities(model.getProp(state,names{i}));
                state = model.setProp(state,names{i},val);
            end
            
            state = setupBCSources(model, state);
            
            %% Update Reaction Coupling variables
            
            names={{'NegativeElectrode'},{'PositiveElectrode'}};
            for i=1:numel(names)
                submodel=model.getSubmodel(names{i});
                val = submodel.updateReactionRate(model.getProp(state,names{i}));
                state = model.setProp(state,names{i},val);
            end
            
            state = setupExchanges(model, state);
            
            %% Update Fluxes variables
            
            names={{'Electrolyte'},{'NegativeElectrode'},{'PositiveElectrode'}};
            for i=1:numel(names)
                submodel=model.getSubmodel(names{i});
                val = submodel.updateQuantities(model.getProp(state,names{i}));
                state = model.setProp(state,names{i},val);
            end
            
            names={{'PositiveCurrentCollector'},{'NegativeCurrentCollector'}};
            for i=1:numel(names)
                submodel=model.getSubmodel(names{i});
                val = submodel.updateChargeCont(model.getProp(state,names{i}));
                state = model.setProp(state,names{i},val);
            end
            
            %% Set up the governing equations
            
            %% Mass and charge conservation for the electorlyte and the electrode
            
            % Accumulation terms for the mass conservation equtions
            cdotLi = struct();
            cdotLi.Electrolyte = (state.Electrolyte.cs{1} - state0.Electrolyte.cs{1})/dt;
            cdotLi.NegativeElectrode = (state.NegativeElectrode.ActiveMaterial.Li - state0.NegativeElectrode.ActiveMaterial.Li)/dt;
            cdotLi.PositiveElectrode = (state.PositiveElectrode.ActiveMaterial.Li - state0.PositiveElectrode.ActiveMaterial.Li)/dt;
            
            names = {'Electrolyte','NegativeElectrode','PositiveElectrode'};
            eqs={};
            
            for i = 1 : numel(names)
                
                submodel = model.getSubmodel({names{i}});

                %% probably only be done on the submodel
                source = model.getProp(state,{names{i},'LiSource'});
                flux = model.getProp(state,{names{i},'LiFlux'});

                %% could use submodel
                div = submodel.operators.Div(flux)./submodel.G.cells.volumes;

                if strcmp(names{i}, 'Electrolyte')
                    cepsdot = submodel.volumeFraction.*cdotLi.(names{i});
                else
                    cepsdot = submodel.ActiveMaterial.volumeFraction.*cdotLi.(names{i});
                end
                %% Li conservation
                eqs{end+1} = div - source + cepsdot;
                
                %¤ charge continutity
                %% should probably be done on the sub model
                eqs{end+1} = model.getProp(state,{names{i},'chargeCont'});
            end
            
            %% charge conservation for the current collectors
            
            names = {'NegativeCurrentCollector','PositiveCurrentCollector'};
            for i = 1 : numel(names)
                eqs{end+1} = model.getProps(state, {names{i}, 'chargeCont'});
            end
            
            %% setup control equation (fixed total current at PositiveCurrentCollector)
            
            src = drivingForces.src(time);
            coupterm = model.getCoupTerm('bc-PositiveCurrentCollector');
            faces = coupterm.couplingfaces;
            bcval = state.PositiveCurrentCollector.E;
            cond_pcc = model.PositiveCurrentCollector.effectiveElectronicConductivity;
            [trans_pcc, cells] = model.PositiveCurrentCollector.operators.harmFaceBC(cond_pcc, faces);
            control = src - sum(trans_pcc.*(bcval - state.PositiveCurrentCollector.phi(cells)));
            
            eqs{end+1} = -control;
            
            
            %% Give type and names to equations and names of the primary variables (for book-keeping)
            
            types={'cell','cell','cell','cell',...
                   'cell','cell','cell','cell','cell'};
            names = {'Electrolyte_Li_massCont', ...
                     'Electrolyte_chargeCont' , ...
                     'NegativeElectrode_Li_massCont'   , ...
                     'NegativeElectrode_e_chargeCont'  , ...
                     'PositiveElectrode_Li_massCont'   , ...
                     'PositiveElectrode_e_chargeCont'  , ...
                     'NegativeCurrentCollector_e_chargeCont', ...
                     'PositiveCurrentCollector_e_chargeCont', ...
                     'control'};
            primaryVars = model.getPrimaryVariables();

            %% setup LinearizedProblem that can be processed by MRST Newton API
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
            
        end
        
        function state = dispatchValues(model, state)
        % Abbreviations used in this function:
        % elyte : Electrolyte
        % ne    : NegativeElectrode
        % pe    : PositiveElectrode
        % ncc   : NegativeCurrentCollector
        % pcc   : PositiveCurrentCollector
            
            T = state.T;
            SOC = state.SOC;

            elyte = model.Electrolyte;
            G = elyte.G;
            T_elyte = T(G.mappings.cellmap);

            ne = model.NegativeElectrode;
            G = ne.G;
            T_ne = T(G.mappings.cellmap);
            SOC_ne = SOC(G.mappings.cellmap);
            
            pe = model.PositiveElectrode;
            G = pe.G;
            T_pe = T(G.mappings.cellmap);
            SOC_pe = SOC(G.mappings.cellmap);
            
            pcc = model.PositiveCurrentCollector;
            G = pcc.G;
            T_pcc = T(G.mappings.cellmap);

            ncc = model.NegativeCurrentCollector;
            G = ncc.G;
            T_ncc = T(G.mappings.cellmap);

            state.Electrolyte.T                      = T_elyte;
            state.NegativeElectrode.T                = T_ne;
            state.NegativeElectrode.ActiveMaterial.T = T_ne;
            state.PositiveElectrode.T                = T_pe;
            state.PositiveElectrode.ActiveMaterial.T = T_pe;
            state.PositiveCurrentCollector.T         = T_pcc;
            state.NegativeCurrentCollector.T         = T_ncc;
            
            state.NegativeElectrode.SOC                = SOC_ne;
            state.NegativeElectrode.ActiveMaterial.SOC = SOC_ne;
            state.PositiveElectrode.SOC                = SOC_pe;
            state.PositiveElectrode.ActiveMaterial.SOC = SOC_pe;
            
        end

        function state = updatePhiElectrolyte(model, state)
        % Abbreviations used in this function 
        % elyte : Electrolyte
        % ne    : NegativeElectrode
        % pe    : PositiveElectrode
            
            
            phi_elyte = state.Electrolyte.phi;
            elyte = model.Electrolyte;

            elyte_cells = zeros(model.G.cells.num, 1);
            elyte_cells(elyte.G.mappings.cellmap) = (1 : elyte.G.cells.num)';

            ne = model.NegativeElectrode;
            phi_elyte_ne = phi_elyte(elyte_cells(ne.G.mappings.cellmap));

            pe = model.PositiveElectrode;
            phi_elyte_pe = phi_elyte(elyte_cells(pe.G.mappings.cellmap));

            state.NegativeElectrode.phiElectrolyte = phi_elyte_ne;
            state.PositiveElectrode.phiElectrolyte = phi_elyte_pe;

        end

        
        function model = setElectrolytePorosity(model)
        % Abbreviations used in this function 
        % elyte : Electrolyte
        % ne    : NegativeElectrode
        % pe    : PositiveElectrode
            
            elyte = model.getSubmodel({'Electrolyte'});
            ne = model.getSubmodel({'NegativeElectrode'});
            pe = model.getSubmodel({'PositiveElectrode'});
            sep = model.getSubmodel({'sep'});

            elyte_cells = zeros(model.G.cells.num, 1);
            elyte_cells(elyte.G.mappings.cellmap) = (1 : elyte.G.cells.num)';

            elyte.volumeFraction = NaN(elyte.G.cells.num, 1);
            elyte.volumeFraction(elyte_cells(ne.G.mappings.cellmap))  = ne.porosity;
            elyte.volumeFraction(elyte_cells(pe.G.mappings.cellmap))  = pe.porosity;
            elyte.volumeFraction(elyte_cells(sep.G.mappings.cellmap)) = sep.porosity;

            model.Electrolyte = elyte;

        end
        
        function initstate = setupInitialState(model)
        % Setup initial state
        %
        % Abbreviations used in this function 
        % elyte : Electrolyte
        % ne    : NegativeElectrode
        % pe    : PositiveElectrode
        % ncc   : NegativeCurrentCollector
        % pcc   : PositiveCurrentCollector
            
            nc = model.G.cells.num;

            SOC = model.SOC;
            T   = model.T;
            
            initstate.T   =  T*ones(nc, 1);
            initstate.SOC =  SOC*ones(nc, 1);

            initstate = model.dispatchValues(initstate);
            
            elyte = model.Electrolyte;
            ne    = model.NegativeElectrode;
            pe    = model.PositiveElectrode;
            ncc   = model.NegativeCurrentCollector;
            pcc   = model.PositiveCurrentCollector;
            ne_am = ne.ActiveMaterial;
            pe_am = pe.ActiveMaterial;

            %% setup initial NegativeElectrode state

            m = (1 ./ (ne_am.theta100 - ne_am.theta0));
            b = -m .* ne_am.theta0;
            theta = (SOC - b) ./ m;
            c = theta .* ne_am.Li.cmax;
            c = c*ones(ne.G.cells.num, 1);

            initstate.NegativeElectrode.ActiveMaterial.Li = c;
            initstate.NegativeElectrode.ActiveMaterial = ne_am.updateQuantities(initstate.NegativeElectrode.ActiveMaterial);

            OCP = initstate.NegativeElectrode.ActiveMaterial.OCP;
            initstate.NegativeElectrode.ActiveMaterial.phi = OCP;

            %% setup initial PositiveElectrode state

            m = (1 ./ (pe_am.theta100 - pe_am.theta0));
            b = -m .* pe_am.theta0;
            theta = (SOC - b) ./ m;
            c = theta .* pe_am.Li.cmax;
            c = c*ones(pe.G.cells.num, 1);

            initstate.PositiveElectrode.ActiveMaterial.Li = c;
            initstate.PositiveElectrode.ActiveMaterial = pe_am.updateQuantities(initstate.PositiveElectrode.ActiveMaterial);

            OCP = initstate.PositiveElectrode.ActiveMaterial.OCP;
            initstate.PositiveElectrode.ActiveMaterial.phi = OCP;

            %% setup initial Electrolyte state

            initstate.Electrolyte.phi = zeros(elyte.G.cells.num, 1);
            cs=cell(2,1);
            initstate.Electrolyte.cs = cs;
            initstate.Electrolyte.cs{1} = 1000*ones(elyte.G.cells.num, 1);

            %% setup initial Current collectors state
            
            OCP = initstate.NegativeElectrode.ActiveMaterial.OCP;
            OCP = OCP(1) .* ones(ncc.G.cells.num, 1);
            initstate.NegativeCurrentCollector.phi = OCP;

            OCP = initstate.PositiveElectrode.ActiveMaterial.OCP;
            OCP = OCP(1) .* ones(pcc.G.cells.num, 1);
            initstate.PositiveCurrentCollector.phi = OCP;
            
            initstate.PositiveCurrentCollector.E = OCP(1);
            
        end
        
        function coupterm = getCoupTerm(model, coupname)
            coupnames = model.couplingnames;
            
            [isok, ind] = ismember(coupname, coupnames);
            assert(isok, 'name of coupling term is not recognized.');
            
            coupterm = model.couplingTerms{ind};
            
        end
        
        function state = initStateAD(model,state)
            adbackend = model.AutoDiffBackend();
            [state.Electrolyte.cs{1},...
             state.Electrolyte.phi,...   
             state.NegativeElectrode.ActiveMaterial.Li,...    
             state.NegativeElectrode.ActiveMaterial.phi,...   
             state.PositiveElectrode.ActiveMaterial.Li,...    
             state.PositiveElectrode.ActiveMaterial.phi,...   
             state.NegativeCurrentCollector.phi,...    
             state.PositiveCurrentCollector.phi,...    
             state.PositiveCurrentCollector.E]=....
                adbackend.initVariablesAD(...
                    state.Electrolyte.cs{1},...
                    state.Electrolyte.phi,...   
                    state.NegativeElectrode.ActiveMaterial.Li,...    
                    state.NegativeElectrode.ActiveMaterial.phi,...   
                    state.PositiveElectrode.ActiveMaterial.Li,...    
                    state.PositiveElectrode.ActiveMaterial.phi,...   
                    state.NegativeCurrentCollector.phi,...    
                    state.PositiveCurrentCollector.phi,...    
                    state.PositiveCurrentCollector.E);       
        end
        
        function p = getPrimaryVariables(model)
            p ={{'Electrolyte','cs'},...
                {'Electrolyte','phi'},...   
                {'NegativeElectrode','ActiveMaterial','Li'},...    
                {'NegativeElectrode','ActiveMaterial','phi'},...   
                {'PositiveElectrode','ActiveMaterial','Li'},...    
                {'PositiveElectrode','ActiveMaterial','phi'},...   
                {'NegativeCurrentCollector','phi'},...    
                {'PositiveCurrentCollector','phi'},...    
                {'PositiveCurrentCollector','E'}
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
            state.Electrolyte.cs{1} =  state.Electrolyte.cs{1} + dx{1};
            report = [];
        end
        
        
        function state = reduceState(model, state, removeContainers)
        % Reduce state to double (do not use property containers)
            state = value(state);
        end

    end

end
