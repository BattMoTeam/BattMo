classdef BatteryModel < CompositeModel

    properties
        fv
    end
    
    methods
        
        function model = BatteryModel(varargin)
            
            model = model@CompositeModel('battery');
            
            G = cartGrid([5, 5]);
            G = computeGeometry(G);
            
            model.G = G;

            nc = G.cells.num;
            cells = (1 : nc)';
            
            sepnx  = 10;
            nenx   = 10;
            penx   = 10;
            ccnenx = 5;
            ccpenx = 5;
            
            nxs = [ccnenx; nenx; sepnx; penx; ccpenx];
            ny = 10;

            xlength = 1e-6*ones(5, 1);
            ylength = 1e-6;
            
            x = xlength./nxs;
            x = rldecode(x, nxs);
            x = [0; cumsum(x)];
            
            y = ylength/ny;
            y = rldecode(y, ny);
            y = [0; cumsum(y)];
            
            G = tensorGrid(x, y);
            G = computeGeometry(G);
            obj.G = G;
            
            obj.componentnames = {'elyte', 'ne', 'pe', 'ccne', 'ccpe'};
            
            %% setup elyte
            nx = sum(nxs); 
            
            submodels = {};
            submodelnames = {};
            
            istart = ccnenx + 1;
            ni = nenx + sepnx + penx;
            cells = pickTensorCells(istart, ni, nx, ny);
            submodels{end + 1} = orgLiPF6('elyte', G, cells);
            
            %% setup ne
            istart = ccnenx + 1;
            cells = pickTensorCells(istart, nenx, nx, ny);
            submodels{end + 1} = graphiteElectrode('ne', G, cells);
            
            %% setup pe
            istart = ccnenx + nenx + sepnx + 1;
            cells = pickTensorCells(istart, penx, nx, ny);
            submodels{end + 1} = nmc111Electrode('pe', G, cells);

            %% setup ccne
            istart = 1;
            cells = pickTensorCells(istart, ccnenx, nx, ny);
            submodels{end + 1} = currentCollector('ccne', G, cells);
            
            %% setup ccpe
            istart = ccnenx + nenx + sepnx + penx + 1;
            cells = pickTensorCells(istart, ccpenx, nx, ny);
            submodels{end + 1}  = currentCollector('ccpe', G, cells);

            %% setup sep 
            istart = ccnenx + nenx + 1;
            cells = pickTensorCells(istart, sepnx, nx, ny);
            submodels{end + 1} = celgard2500('sep', G, cells);

            submodelnames = cellfun(@(m) m.getModelName(), submodels, 'uniformoutput', false);
            
            model.SubModels = submodels;
            model.SubModelNames = submodelnames;
            
            model.isnamespaceroot = true;
            model = model.initiateCompositeModel();
            
        end

        function state = initializeState(model, state)
            
            state = model.validateState(state);
            
            % initialize each of the submodels
            for i = 1 : model.nSubModels
                state = model.SubModels{i}.initializeState(state);
            end
            
            ccne = model.getSubModel('ccne');
            ne = model.getSubModel('ne');
            
            OCP_ne = ne.getProp(state, 'OCP');
            nc = ccne.G.cells.num;
            OCP_ccne = OCP_ne(1)*ones(nc, 1);
            
            state = ccne.setProp(state, 'OCP', OCP_ccne);
            state = ccne.setProp(state, 'phi', OCP_ccne);
            
            ccpe = model.getSubModel('ccpe');
            pe = model.getSubModel('pe');
            
            OCP_pe = pe.getProp(state, 'OCP');
            nc = ccpe.G.cells.num;
            OCP_ccpe = OCP_pe(1)*ones(nc, 1);
            
            state = ccpe.setProp(state, 'OCP', OCP_ccpe);
            state = ccpe.setProp(state, 'phi', OCP_ccpe);
            
            state = ccpe.setProp(state, 'E', 0);
            
        end

        function model = setupFV(model, state)
            model.fv = fv2d(model, state);
        end
        
        function [namespaces, names] = getModelPrimaryVarNames(model)
           [namespaces, names] = getModelPrimaryVarNames@CompositeModel(model);
           namespaces = {namespaces{:}, 'ccpe'};
           names = {names{:}, 'E'};
        end
        
        function [namespaces, names] = getModelVarNames(model)
            [namespaces1, names1] = getModelVarNames@CompositeModel(model);
            names2 = {'T', 'SOC'};
            namespaces2 = model.assignCurrentNameSpace(names2);
            
            namespaces = horzcat(namespaces1, namespaces2);
            names = horzcat(names1, names2);
        end
        
        function [y, yp] = dynamicPreprocess(model, state)

            %% Initialize the state vector
            
            y0 = [];
            pvarnames = model.getPrimaryVarNames();
            for ind = 1 : numel(pvarnames);
                y0 = [y0; model.getProp(state, pvarnames{ind})];
            end
            
            yp0  = zeros(length(y0), 1);
            
        end

        
        function dynamicBuildSOE(obj, t, y, yp, varargin)

            opt = struct('useAD', false);
            opt = merge_options(opt, varargin{:});
            useAD = opt.useAD;

            fv = obj.fv;

            if useAD
                adbackend = obj.AutoDiffBackend();
                [y, yp] = adbackend.initVariablesAD(y, yp);
            end

            % Mapping of variables  
            
            % elyte variables
            obj.elyte.sp.Li.ceps = y(fv.getSlot('elyte_Li'));
            obj.elyte.phi = y(fv.getSlot('elyte_phi'));
            % ne variables
            obj.ne.am.Li.cseps = y(fv.getSlot('ne_Li'));
            obj.ne.am.phi = y(fv.getSlot('ne_phi'));
            % pe variables
            obj.pe.am.Li.cseps = y(fv.getSlot('pe_Li'));
            obj.pe.am.phi = y(fv.getSlot('pe_phi'));
            % ccne variables
            obj.ccne.am.phi = y(fv.getSlot('ccne_phi'));
            % ccpe variables
            obj.ccpe.am.phi = y(fv.getSlot('ccpe_phi'));
            % voltage closure variable
            obj.ccpe.E = y(fv.getSlot('E'));
            
            % variables for time derivatives
            obj.elyte.sp.Li.cepsdot = yp(fv.getSlot('elyte_Li'));
            obj.ne.am.Li.csepsdot   = yp(fv.getSlot('ne_Li'));
            obj.pe.am.Li.csepsdot   = yp(fv.getSlot('pe_Li'));

            obj.elyte.sp.Li.c  = obj.elyte.sp.Li.ceps ./ obj.elyte.eps;
            obj.elyte.sp.PF6.c = obj.elyte.sp.Li.c;

            obj.ne.am.Li.cs = obj.ne.am.Li.cseps ./ obj.ne.am.eps;
            obj.pe.am.Li.cs = obj.pe.am.Li.cseps ./ obj.pe.am.eps;

            %% Update electrolyte physicochemical and transport properties
            obj.elyte.update()
            obj.elyte.kappaeff = obj.elyte.kappa .* obj.elyte.eps .^1.5;
            obj.elyte.sp.Li.Deff = obj.elyte.sp.Li.D .* obj.elyte.eps .^1.5;

            obj.ne.am.update()
            obj.ne.am.Li.Deff = obj.ne.am.Li.D .* obj.ne.am.eps.^1.5;
            obj.ne.sigmaeff = obj.ne.am.sigma .* obj.ne.am.eps.^1.5;

            obj.pe.am.update()
            obj.pe.am.Li.Deff = obj.pe.am.Li.D .* obj.pe.am.eps.^1.5;
            obj.pe.sigmaeff = obj.pe.am.sigma .* obj.pe.am.eps.^1.5;

            obj.ccne.sigmaeff = obj.ccne.am.sigma .* obj.ccne.am.eps.^1.5;
            obj.ccpe.sigmaeff = obj.ccpe.am.sigma .* obj.ccpe.am.eps.^1.5;

            %% Ionic current density
            % Ionic current density in the liquid
            %   Ionic current density due to the chemical potential gradient
            N = obj.elyte.N;
            ncomp = obj.elyte.ncomp;

            jchems = cell(ncomp, 1);
            for i = 1 : ncomp
                jchems{i} = zeros(N + 1, 1);
                coeff =  obj.elyte.kappaeff .* obj.elyte.ion.tvec{i} .* obj.elyte.ion.dmudc{i} ./ ...
                                     (obj.elyte.ion.zvec{i}.*obj.con.F);
                jchems{i} = obj.elyte.operators.harmFace(coeff).* obj.elyte.operators.Grad(obj.elyte.ion.cvec{i});

            end
            obj.elyte.jchem = jchems{1};
            for i = 2 : ncomp
                obj.elyte.jchem = obj.elyte.jchem + jchems{i};
            end
            %   Ionic current density due to the electrochemical potential gradient
            obj.elyte.j = obj.elyte.operators.harmFace(obj.elyte.kappaeff).*(-1).*obj.elyte.operators.Grad(obj.elyte.phi) - obj.elyte.jchem;

            %% Electric current density
            % Active material NE

            ccne = obj.ccne;
            ne = obj.ne;

            obj.ne.j = ne.operators.harmFace(ne.sigmaeff) .* (-1) .* ne.operators.Grad(ne.am.phi);
            obj.ccne.j = ccne.operators.harmFace(ccne.sigmaeff) .* (-1) .* ccne.operators.Grad(ccne.am.phi);

            % Add current transfers between ccne collector and ne material. They correspond to flux continuity
            coupterm = obj.getCoupTerm('ccne-ne');
            face_ccne = coupterm.couplingfaces(:, 1);
            face_ne = coupterm.couplingfaces(:, 2);
            [tne, bccell_ne] = ne.operators.harmFaceBC(ne.sigmaeff, face_ne);
            [tccne, bccell_ccne] = ccne.operators.harmFaceBC(ccne.sigmaeff, face_ccne);
            bcphi_ne = ne.am.phi(bccell_ne);
            bcphi_ccne = ccne.am.phi(bccell_ccne);

            obj.ne.j_bcsource = ne.am.phi*0.0; %NB hack to initialize zero ad
            obj.ccne.j_bcsource = ccne.am.phi*0.0; %NB hack to initialize zero ad

            trans = 1./(1./tne + 1./tccne);
            crosscurrent = trans.*(bcphi_ccne - bcphi_ne);
            obj.ne.j_bcsource(bccell_ne) = crosscurrent;
            obj.ccne.j_bcsource(bccell_ccne) = -crosscurrent;

            % We impose the boundary condition at chosen boundary cells of the ne current collector
            coupterm = obj.getCoupTerm('bc-ccne');
            faces = coupterm.couplingfaces;
            bcval = zeros(numel(faces), 1);
            [tccne, cells] = obj.ccne.operators.harmFaceBC(obj.ccne.sigmaeff, faces);
            obj.ccne.j_bcsource(cells) = obj.ccne.j_bcsource(cells) + tccne.*(bcval - obj.ccne.am.phi(cells));

            % Active material PE and current collector

            ccpe = obj.ccpe;
            pe = obj.pe;

            obj.pe.j =  pe.operators.harmFace(pe.sigmaeff) .* (-1) .* pe.operators.Grad(pe.am.phi);
            obj.ccpe.j =  ccpe.operators.harmFace(ccpe.sigmaeff) .* (-1) .* ccpe.operators.Grad(ccpe.am.phi);

            % Add current transfers between ccpe collector and pe material. They correspond to flux continuity
            coupterm = obj.getCoupTerm('ccpe-pe');
            face_ccpe = coupterm.couplingfaces(:, 1);
            face_pe = coupterm.couplingfaces(:, 2);
            [tpe, bccell_pe] = pe.operators.harmFaceBC(pe.sigmaeff, face_pe);
            [tccpe, bccell_ccpe] = ccpe.operators.harmFaceBC(ccpe.sigmaeff, face_ccpe);
            bcphi_pe = pe.am.phi(bccell_pe);
            bcphi_ccpe = ccpe.am.phi(bccell_ccpe);

            obj.pe.j_bcsource   = pe.am.phi*0.0; %NB hack to initialize zero ad
            obj.ccpe.j_bcsource = ccpe.am.phi*0.0; %NB hack to initialize zero ad

            trans = 1./(1./tpe + 1./tccpe);
            crosscurrent = trans.*(bcphi_ccpe - bcphi_pe);
            obj.pe.j_bcsource(bccell_pe) = crosscurrent;
            obj.ccpe.j_bcsource(bccell_ccpe) = -crosscurrent;

            % We impose the boundary condition at chosen boundary cells of the anode current collector
            coupterm = obj.getCoupTerm('bc-ccpe');
            faces = coupterm.couplingfaces;
            bcval = obj.ccpe.E;
            [tccpe, cells] = obj.ccpe.operators.harmFaceBC(obj.ccpe.sigmaeff, faces);
            obj.ccpe.j_bcsource(cells) = obj.ccpe.j_bcsource(cells) + tccpe.*(bcval - obj.ccpe.am.phi(cells));

            %% Cell voltage
            obj.ccne.E = 0;
            obj.U = obj.ccpe.E - obj.ccne.E;
            
            if useAD
                adsample = getSampleAD(y, yp);
                adbackend = obj.AutoDiffBackend;
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Source terms for continuity equations                    %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Inititalize vectors
            obj.elyte.sp.Li.source = zeros(obj.elyte.N, 1);
            obj.ne.am.Li.source    = zeros(obj.ne.N, 1);
            obj.pe.am.Li.source    = zeros(obj.pe.N, 1);
            obj.ne.am.e.source     = zeros(obj.ne.N, 1);
            obj.pe.am.e.source     = zeros(obj.pe.N, 1);

            if useAD
                obj.elyte.sp.Li.source = adbackend.convertToAD(obj.elyte.sp.Li.source, adsample);
                obj.ne.am.Li.source    = adbackend.convertToAD(obj.ne.am.Li.source, adsample);
                obj.pe.am.Li.source    = adbackend.convertToAD(obj.pe.am.Li.source, adsample);
                obj.ne.am.e.source     = adbackend.convertToAD(obj.ne.am.e.source, adsample);
                obj.pe.am.e.source     = adbackend.convertToAD(obj.pe.am.e.source, adsample);
            end

            %%%%% Set up chemical source terms %%%%%%%%%%%%%%%%%k
            
            %%%%% NE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            coupterm = obj.getCoupTerm('ne-elyte');
            necells = coupterm.couplingcells(:, 1);
            elytecells = coupterm.couplingcells(:, 2);

            % calculate rection rate
            obj.ne.reactBV(obj.elyte.phi(elytecells));

            % Electrolyte NE Li+ source
            obj.elyte.sp.Li.source(elytecells) = +1 .* obj.ne.R;
            
            % Active Material NE Li0 source
            obj.ne.am.Li.source(necells) = -1 .* obj.ne.R;
            
            % Active Material NE current source
            obj.ne.am.e.source(necells) = +1 .* obj.ne.R;

            %%%%% PE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Electrolyte PE Li+ source
            coupterm = obj.getCoupTerm('pe-elyte');
            pecells = coupterm.couplingcells(:, 1);
            elytecells = coupterm.couplingcells(:, 2);
            
            % calculate rection rate
            obj.pe.reactBV(obj.elyte.phi(elytecells));

            % Electrolyte PE Li+ source
            obj.elyte.sp.Li.source(elytecells) = -1 .* obj.pe.R;

            % Active Material PE Li0 source
            obj.pe.am.Li.source(pecells) = +1 .* obj.pe.R;

            % Active Material PE current source
            obj.pe.am.e.source(pecells) = -1 .* obj.pe.R;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Diffusion Flux                                           %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Divergence of diffusion mass flux
            % Electrolyte Li+ Diffusion

            x = obj.elyte.sp.Li.c;
            flux = - obj.elyte.sp.Li.Trans.*obj.elyte.operators.Grad(x);
            obj.elyte.sp.Li.divDiff = obj.elyte.operators.Div(flux)./obj.elyte.Grid.cells.volumes;

            x= obj.ne.am.Li.cs;
            flux = - obj.ne.am.Li.Trans.*obj.ne.operators.Grad(x);
            obj.ne.am.Li.divDiff =  obj.ne.operators.Div(flux)./obj.ne.Grid.cells.volumes;

            x = obj.pe.am.Li.cs;
            flux = - obj.pe.am.Li.Trans.*obj.pe.operators.Grad(x);
            obj.pe.am.Li.divDiff = obj.pe.operators.Div(flux)./obj.pe.Grid.cells.volumes;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Migration Flux                                           %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Divergence of the migration mass flux
            %   Electrolyte Li+ Migration
            flux = obj.elyte.sp.Li.t ./ (obj.elyte.sp.Li.z .* obj.con.F) .* obj.elyte.j;

            obj.elyte.sp.Li.divMig = obj.elyte.operators.Div(flux)./obj.elyte.Grid.cells.volumes;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% System of Equations                                      %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %% Liquid electrolyte dissolved ionic species mass continuity %
            %   Electrolyte Li+ Mass Continuity
            obj.elyte.sp.Li.massCont = ( - obj.elyte.sp.Li.divDiff - obj.elyte.sp.Li.divMig + obj.elyte.sp.Li.source ...
                                         - obj.elyte.sp.Li.cepsdot);

            %% Liquid electrolyte charge continuity %%%%%%%%%%%%%%%%%%%%%%%

            obj.elyte.chargeCont = -(obj.elyte.operators.Div( obj.elyte.j)./obj.elyte.Grid.cells.volumes) ./ obj.con.F+ ...
                obj.elyte.sp.Li.source .* obj.elyte.sp.Li.z;

            %% Active material mass continuity %%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.ne.am.Li.massCont = (-obj.ne.am.Li.divDiff + obj.ne.am.Li.source - obj.ne.am.Li.csepsdot);
            obj.pe.am.Li.massCont = (-obj.pe.am.Li.divDiff + obj.pe.am.Li.source - obj.pe.am.Li.csepsdot);

            %% Active material charge continuity %%%%%%%%%%%%%%%%%%%%%%%%%%%

            obj.ne.am.e.chargeCont = (obj.ne.operators.Div(obj.ne.j) - obj.ne.j_bcsource)./ ...
                obj.ne.Grid.cells.volumes./obj.con.F - obj.ne.am.e.source;
            obj.pe.am.e.chargeCont = (obj.pe.operators.Div(obj.pe.j) - obj.pe.j_bcsource)./ ...
                obj.pe.Grid.cells.volumes./obj.con.F - obj.pe.am.e.source;

            obj.ccne.am.e.chargeCont = (obj.ccne.operators.Div(obj.ccne.j) - obj.ccne.j_bcsource)./ ...
                obj.ccne.Grid.cells.volumes./obj.con.F;
            obj.ccpe.am.e.chargeCont = (obj.ccpe.operators.Div(obj.ccpe.j) - obj.ccpe.j_bcsource)./ ...
                obj.ccpe.Grid.cells.volumes./obj.con.F;

            %% control equation
            src = currentSource(t, fv.tUp, fv.tf, obj.J);
            coupterm = obj.getCoupTerm('bc-ccpe');
            faces = coupterm.couplingfaces;
            bcval = obj.ccpe.E;
            [tccpe, cells] = obj.ccpe.operators.harmFaceBC(obj.ccpe.sigmaeff, faces);
            control = src - sum(tccpe.*(bcval - obj.ccpe.am.phi(cells)));

            %% State vector %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.soe = vertcat(obj.elyte.sp.Li.massCont, ...
                              obj.elyte.chargeCont    , ...
                              obj.ne.am.Li.massCont   , ...
                              obj.ne.am.e.chargeCont  , ...
                              obj.pe.am.Li.massCont   , ...
                              obj.pe.am.e.chargeCont  , ...
                              obj.ccne.am.e.chargeCont, ...
                              obj.ccpe.am.e.chargeCont, ...
                              control);

        end

    end
    
    
end
