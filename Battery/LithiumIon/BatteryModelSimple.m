classdef BatteryModelSimple < PhysicalModel %< CompositeModel

    properties

        con = physicalConstants();

        % Coupling structures
        couplingnames
        couplingTerms

        % fv2d structure (will disappear when we switch to own Newton solver)
        fv

        % Temperature and SOC
        % for the moment here, for convenience. Will be moved
        T
        SOC

        % Input current
        J
        % Voltage cut
        Ucut

        % submodels
        elyte
        ne
        pe
        ccpe
        ccne
        sep
        
    end

    methods

        function model = BatteryModelSimple(varargin)

            model = model@PhysicalModel([]);
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks',false);
            names = {'T', 'SOC'};
            %model.names = names;
            %model = model.setupVarDims();

            model = model.setupBatteryComponents();
            model = model.setElytePorosity();
            %% setup couplings
            coupTerms = {};

            % coupling term 'ne-cc'
            coupTerms{end + 1} = setupNeElyteCoupTerm(model);
            coupTerms{end + 1} = setupPeElyteCoupTerm(model);
            coupTerms{end + 1} = setupCcneNeCoupTerm(model);
            coupTerms{end + 1} = setupCcpePeCoupTerm(model);
            coupTerms{end + 1} = setupCcneBcCoupTerm(model);
            coupTerms{end + 1} = setupCcpeBcCoupTerm(model);

            model.couplingTerms = coupTerms;
            model.couplingnames = cellfun(@(x) x.name, coupTerms, 'uniformoutput', false);

            model.SOC = 0.5;
            model.T = 298.15;

            model.J = 0.1;
            model.Ucut = 2;
            
        end
        
        
        function model = setupBatteryComponents(model)
            
            fac = 1;
            sepnx  = 30*fac;
            nenx   = 30*fac;
            penx   = 30*fac;
            ccnenx = 20*fac;
            ccpenx = 20*fac;

            nxs = [ccnenx; nenx; sepnx; penx; ccpenx];
            ny = 10*fac;

            xlength = 1e-6*[10; 100; 50; 80; 10];
            ylength = 1e-2;

            x = xlength./nxs;
            x = rldecode(x, nxs);
            x = [0; cumsum(x)];

            y = ylength/ny;
            y = rldecode(y, ny);
            y = [0; cumsum(y)];

            G = tensorGrid(x, y);
            G = computeGeometry(G);
            model.G = G;

            %% setup elyte
            nx = sum(nxs);



            istart = ccnenx + 1;
            ni = nenx + sepnx + penx;
            cells = pickTensorCells(istart, ni, nx, ny);
            %submodels{end + 1} = orgLiPF6('elyte', G, cells);
            model.elyte = orgLiPF6('elyte', G, cells);
            %% setup ne
            istart = ccnenx + 1;
            cells = pickTensorCells(istart, nenx, nx, ny);
            %submodels{end + 1} = graphiteElectrode('ne', G, cells);
            model.ne = graphiteElectrode('ne', G, cells);
            %% setup pe
            istart = ccnenx + nenx + sepnx + 1;
            cells = pickTensorCells(istart, penx, nx, ny);
            %submodels{end + 1} = nmc111Electrode('pe', G, cells);
            model.pe = nmc111Electrode('pe', G, cells);
            %% setup ccne
            istart = 1;
            cells = pickTensorCells(istart, ccnenx, nx, ny);
            %submodels{end + 1} = currentCollector('ccne', G, cells);
            model.ccne = currentCollector('ccne', G, cells);
            %% setup ccpe
            istart = ccnenx + nenx + sepnx + penx + 1;
            cells = pickTensorCells(istart, ccpenx, nx, ny);
            %submodels{end + 1}  = currentCollector('ccpe', G, cells);
            model.ccpe = currentCollector('ccpe', G, cells);
            %% setup sep
            istart = ccnenx + nenx + 1;
            cells = pickTensorCells(istart, sepnx, nx, ny);
            %submodels{end + 1} = celgard2500('sep', G, cells);
            model.sep = celgard2500('sep', G, cells);
            %model.SubModels = submodels;
            
        end
        
        function initstate = icp2d(model)
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

            %elyte = model.getAssocModel('elyte');
            %ne = model.getAssocModel('ne');
            %pe = model.getAssocModel('pe');
            
            elyte = model.getSubmodel({'elyte'});
            ne = model.getSubmodel({'ne'});
            pe = model.getSubmodel({'pe'})
            sep = model.getSubmodel({'sep'});

            elytecells = zeros(model.G.cells.num, 1);
            elytecells(elyte.G.mappings.cellmap) = (1 : elyte.G.cells.num)';

            elyte.eps = NaN(elyte.G.cells.num, 1);
            elyte.eps(elytecells(ne.G.mappings.cellmap)) = ne.void;
            elyte.eps(elytecells(pe.G.mappings.cellmap)) = pe.void;
            elyte.eps(elytecells(sep.G.mappings.cellmap)) = sep.void;
            model.elyte =elyte;
            %model = model.setSubModel('elyte', elyte);

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
         
        function [problem, state] = getEquations(model, state0, state,dt, drivingForces, varargin)
            time = state0.time+dt;
            state=model.initStateAD(state);
            
            nc = model.G.cells.num;
            
            % setup temperature and SOC here
            %% for now this is kept constant?
            state.T =  model.T*ones(nc, 1);
            state.SOC =  model.SOC*ones(nc, 1);
            
            % variables for time derivatives
            cdotLi=struct();
            cdotLi.elyte = (state.elyte.cs{1} - state0.elyte.cs{1})/dt;
            cdotLi.ne    = (state.ne.am.Li - state0.ne.am.Li)/dt;
            cdotLi.pe   = (state.pe.am.Li - state0.pe.am.Li)/dt;
            
            state = model.dispatchValues(state);
            state = model.updatePhiElyte(state);
            %% first update level 2
            names={{'pe','am'},{'pe','am'}};
            for i=1:numel(names)
                submodel = model.getSubmodel(names{i});
                val = submodel.updateQuantities(model.getProp(state,names{i}));
                state = model.setProp(state,names{i},val);
            end
            
            
            state = setupBCSources(model, state);
            
            names={{'ne'},{'pe'}};
            for i=1:numel(names)
                submodel=model.getSubmodel(names{i});
                val = submodel.updateReactionRate(model.getProp(state,names{i}));
                state = model.setProp(state,names{i},val);
            end
            state = setupExchanges(model, state);
            
            %%update level 1
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
            
            %% set equations
            names={'elyte','ne','pe'};
            eqs={};
            for i=1:numel(names)
                submodel=model.getSubmodel({names{i}});
                %% probably only be done on the submodel
                source = model.getProp(state,{names{i},'LiSource'});
                flux = model.getProp(state,{names{i},'LiFlux'});
                %% could use submodel
                div =  submodel.operators.Div(flux)./submodel.G.cells.volumes;
                %% HAC
                if(strcmp(names{i},'elyte'))
                    cepsdot = submodel.eps.*cdotLi.(names{i});
                else
                    cepsdot = submodel.am.eps.*cdotLi.(names{i});
                end
                %% Li conservation
                eqs{end+1} = -div + source - cepsdot;
                % charge continutity
                %% should probably be done on the sub model
                eqs{end+1} = model.getProp(state,{names{i},'chargeCont'});
            end
            names={'ccne','ccpe'};
            for i=1:numel(names)
                eqs{end+1} = model.getProps(state,{names{i},'chargeCont'});
            end
            
            %src = currentSource(t, fv.tUp, fv.tf, model.J);
            src = drivingForces.src(time);%%(t, fv.tUp, fv.tf, model.J);
            coupterm = model.getCoupTerm('bc-ccpe');
            faces = coupterm.couplingfaces;
            bcval = state.ccpe.E;
            ccpe_sigmaeff = model.ccpe.sigmaeff;
            [tccpe, cells] = model.ccpe.operators.harmFaceBC(ccpe_sigmaeff, faces);
            control = src - sum(tccpe.*(bcval - state.ccpe.phi(cells)));
            
            %% Governing equations
            
            eqs{end+1} = control;
            
            
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
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
            %state.cdotLi=cdotLi;
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
        % Reduce state to doubles, and optionally remove the property
        % containers to reduce storage space only p = model.getPrimaryVariables();
        % should be need
            state = value(state);
            
        end
        
        function [convergence, values, names] = checkConvergence(model, problem, varargin)
            
            [values, tolerances, names] = getConvergenceValues(model, problem, varargin{:});
            convergence = values < tolerances;
            %if model.verbose
            %    fprintf('Iteration %i ',problem.iterationNo);
            %    fprintf(' residual ');
            %    fprintf(' %d ',values);
            %    fprintf('\n');
            %    disp(values)
            %end
        end
        
        function coupterm = getCoupTerm(model, coupname)
            coupnames = model.couplingnames;
            
            [isok, ind] = ismember(coupname, coupnames);
            assert(isok, 'name of coupling term is not recognized.');
            
            coupterm = model.couplingTerms{ind};
            
        end
        
        function coupTerm = setupNeElyteCoupTerm(model)
            
            ne = model.ne;
            elyte = model.elyte;
            
            Gne = ne.G;
            Gelyte = elyte.G;
            
            % parent Grid
            G = Gne.mappings.parentGrid;
            
            % All the cells from ne are coupled with elyte
            cells1 = (1 : Gne.cells.num)';
            pcells = Gne.mappings.cellmap(cells1);
            
            mapping = zeros(G.cells.num, 1);
            mapping(Gelyte.mappings.cellmap) = (1 : Gelyte.cells.num)';
            cells2 = mapping(pcells);
            
            compnames = {'ne', 'elyte'};
            coupTerm = couplingTerm('ne-elyte', compnames);
            coupTerm.couplingcells =  [cells1, cells2];
            coupTerm.couplingfaces = []; % no coupling throug faces. We set it as empty
            
        end
        
        function coupTerm = setupPeElyteCoupTerm(model)
            
            pe = model.pe;
            elyte = model.elyte;
            
            Gpe = pe.G;
            Gelyte = elyte.G;
            
            % parent Grid
            G = Gpe.mappings.parentGrid;
            
            % All the cells from pe are coupled with elyte
            cells1 = (1 : Gpe.cells.num)';
            pcells = Gpe.mappings.cellmap(cells1);
            
            mapping = zeros(G.cells.num, 1);
            mapping(Gelyte.mappings.cellmap) = (1 : Gelyte.cells.num)';
            cells2 = mapping(pcells);
            
            compnames = {'pe', 'elyte'};
            coupTerm = couplingTerm('pe-elyte', compnames);
            coupTerm.couplingcells = [cells1, cells2];
            coupTerm.couplingfaces = []; % no coupling between faces
            
        end

        function coupTerm = setupCcneNeCoupTerm(model)

            ne = model.ne;
            ccne = model.ccne;

            Gne = ne.G;
            Gccne = ccne.G;

            G = Gne.mappings.parentGrid;

            netbls = setupSimpleTables(Gne);
            ccnetbls = setupSimpleTables(Gccne);
            tbls = setupSimpleTables(G);            
            
            necelltbl = netbls.celltbl;
            necelltbl = necelltbl.addInd('globcells', Gne.mappings.cellmap);
            nefacetbl = netbls.facetbl;
            nefacetbl = nefacetbl.addInd('globfaces', Gne.mappings.facemap);

            necellfacetbl = netbls.cellfacetbl;
            necellfacetbl = crossIndexArray(necellfacetbl, necelltbl, {'cells'});
            necellfacetbl = crossIndexArray(necellfacetbl, nefacetbl, {'faces'});
            
            ccnecelltbl = ccnetbls.celltbl;
            ccnecelltbl = ccnecelltbl.addInd('globcells', Gccne.mappings.cellmap);
            ccnefacetbl = ccnetbls.facetbl;
            ccnefacetbl = ccnefacetbl.addInd('globfaces', Gccne.mappings.facemap);
            
            ccnecellfacetbl = ccnetbls.cellfacetbl;
            ccnecellfacetbl = crossIndexArray(ccnecellfacetbl, ccnecelltbl, {'cells'});
            ccnecellfacetbl = crossIndexArray(ccnecellfacetbl, ccnefacetbl, {'faces'});
            
            gen = CrossIndexArrayGenerator();
            gen.tbl1 = necellfacetbl;
            gen.tbl2 = ccnecellfacetbl;
            gen.replacefds1 = {{'cells', 'necells'}, {'faces', 'nefaces'}, {'globcells', 'neglobcells'}};
            gen.replacefds2 = {{'cells', 'ccnecells'}, {'faces', 'ccnefaces'}, {'globcells', 'ccneglobcells'}};
            gen.mergefds = {'globfaces'};
            
            cell12facetbl = gen.eval();

            ccnefaces = cell12facetbl.get('ccnefaces');
            nefaces = cell12facetbl.get('nefaces');
            ccnecells = cell12facetbl.get('ccnecells');
            necells = cell12facetbl.get('necells');            
            
            compnames = {'ccne', 'ne'};
            coupTerm = couplingTerm('ccne-ne', compnames);
            coupTerm.couplingfaces =  [ccnefaces, nefaces];
            coupTerm.couplingcells = [ccnecells, necells];

        end

        function coupTerm = setupCcpePeCoupTerm(model)

            pe = model.pe;
            ccpe = model.ccpe;

            Gpe = pe.G;
            Gccpe = ccpe.G;

            G = Gpe.mappings.parentGrid;

            petbls = setupSimpleTables(Gpe);
            ccpetbls = setupSimpleTables(Gccpe);
            tbls = setupSimpleTables(G);            
            
            pecelltbl = petbls.celltbl;
            pecelltbl = pecelltbl.addInd('globcells', Gpe.mappings.cellmap);
            pefacetbl = petbls.facetbl;
            pefacetbl = pefacetbl.addInd('globfaces', Gpe.mappings.facemap);

            pecellfacetbl = petbls.cellfacetbl;
            pecellfacetbl = crossIndexArray(pecellfacetbl, pecelltbl, {'cells'});
            pecellfacetbl = crossIndexArray(pecellfacetbl, pefacetbl, {'faces'});
            
            ccpecelltbl = ccpetbls.celltbl;
            ccpecelltbl = ccpecelltbl.addInd('globcells', Gccpe.mappings.cellmap);
            ccpefacetbl = ccpetbls.facetbl;
            ccpefacetbl = ccpefacetbl.addInd('globfaces', Gccpe.mappings.facemap);
            
            ccpecellfacetbl = ccpetbls.cellfacetbl;
            ccpecellfacetbl = crossIndexArray(ccpecellfacetbl, ccpecelltbl, {'cells'});
            ccpecellfacetbl = crossIndexArray(ccpecellfacetbl, ccpefacetbl, {'faces'});
            
            gen = CrossIndexArrayGenerator();
            gen.tbl1 = pecellfacetbl;
            gen.tbl2 = ccpecellfacetbl;
            gen.replacefds1 = {{'cells', 'pecells'}, {'faces', 'pefaces'}, {'globcells', 'peglobcells'}};
            gen.replacefds2 = {{'cells', 'ccpecells'}, {'faces', 'ccpefaces'}, {'globcells', 'ccpeglobcells'}};
            gen.mergefds = {'globfaces'};
            
            cell12facetbl = gen.eval();

            ccpefaces = cell12facetbl.get('ccpefaces');
            pefaces = cell12facetbl.get('pefaces');
            ccpecells = cell12facetbl.get('ccpecells');
            pecells = cell12facetbl.get('pecells');            
            
            compnames = {'ccpe', 'pe'};
            coupTerm = couplingTerm('ccpe-pe', compnames);
            coupTerm.couplingfaces =  [ccpefaces, pefaces];
            coupTerm.couplingcells = [ccpecells, pecells];
            
        end

        function coupTerm = setupCcneBcCoupTerm(model)

            ccne = model.ccne;
            G = ccne.G;

            % We pick up the faces at the top of Cccne
            yf = G.faces.centroids(:, 2);
            myf = max(yf);
            faces = find(abs(yf-myf) < eps*1000);
            cells = sum(G.faces.neighbors(faces, :), 2);

            compnames = {'ccne'};
            coupTerm = couplingTerm('bc-ccne', compnames);
            coupTerm.couplingfaces = faces;
            coupTerm.couplingcells = cells;

        end

        function coupTerm = setupCcpeBcCoupTerm(model)

            ccpe = model.ccpe;
            G = ccpe.G;

            % We pick up the faces at the top of Cccpe
            yf = G.faces.centroids(:, 2);
            myf = max(yf);
            faces = find(abs(yf-myf) < eps*1000q);
            cells = sum(G.faces.neighbors(faces, :), 2);

            compnames = {'ccpe'};
            coupTerm = couplingTerm('bc-ccpe', compnames);
            coupTerm.couplingfaces = faces;
            coupTerm.couplingcells = cells;

        end

    end

end
