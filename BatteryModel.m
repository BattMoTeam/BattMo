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

            model.hasparent = false;

            ccpe = model.getSubModel('ccpe');
            ccpe.pnames = {ccpe.pnames{:}, 'E'};
            model = model.setSubModel(ccpe, 'ccpe');
            
            names = {'T', 'SOC'};
            varnames = model.assignCurrentNameSpace(names);
            
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
            
            OCP_ne = ne.getProp(state, {'graphite', 'OCP'});
            nc = ccne.G.cells.num;
            OCP_ccne = OCP_ne(1)*ones(nc, 1);
            
            state = ccne.setProp(state, 'OCP', OCP_ccne);
            state = ccne.setProp(state, 'phi', OCP_ccne);
            
            ccpe = model.getSubModel('ccpe');
            pe = model.getSubModel('pe');
            
            OCP_pe = pe.getProp(state, {'nmc111', 'OCP'});
            nc = ccpe.G.cells.num;
            OCP_ccpe = OCP_pe(1)*ones(nc, 1);
            
            state = ccpe.setProp(state, 'OCP', OCP_ccpe);
            state = ccpe.setProp(state, 'phi', OCP_ccpe);
            
            state = ccpe.setProp(state, 'E', OCP_ccpe(1));
            
        end

        function model = setupFV(model, state)
            model.fv = fv2d(model, state);
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

        
        function dynamicBuildSOE(model, t, y, yp, varargin)

            opt = struct('useAD', false);
            opt = merge_options(opt, varargin{:});
            useAD = opt.useAD;

            fv = model.fv;

            if useAD
                adbackend = obj.AutoDiffBackend();
                [y, yp] = adbackend.initVariablesAD(y, yp);
            end

            % Mapping of variables  
            state = [];
            state = model.validateState(state);
            
            % short cut
            fun = @(state, name) (model.setProp(state, name, y(fv.getSlot(name))));
            
            % elyte variable
            state = fun(state, 'elyte_c_Li');
            state = fun(state, 'elyte_phi');
            % ne variables
            state = fun(state, 'ne_Li');
            state = fun(state, 'ne_phi');
            % pe variables
            state = fun(state, 'pe_Li');
            state = fun(state, 'pe_phi');
            % ccne variables
            state = fun(state, 'ccne_phi');
            % ccpe variables
            state = fun(state, 'ccpe_phi');
            % voltage closure variable
            state = fun(state, 'ccpe_E');
            
            % variables for time derivatives
            obj.elyte.sp.Li.cepsdot = yp(fv.getSlot('elyte_Li'));
            obj.ne.am.Li.csepsdot   = yp(fv.getSlot('ne_Li'));
            obj.pe.am.Li.csepsdot   = yp(fv.getSlot('pe_Li'));

            % elyte variable
            elyte_c_Li = model.getProp(state, 'elyte_c_Li');
            elyte_phi  = model.getProp(state, 'elyte_phi');
            ne_Li      = model.getProp(state, 'ne_Li');
            ne_phi     = model.getProp(state, 'ne_phi');
            pe_Li      = model.getProp(state, 'pe_Li');
            pe_phi     = model.getProp(state, 'pe_phi');
            ccne_phi   = model.getProp(state, 'ccne_phi');
            ccpe_phi   = model.getProp(state, 'ccpe_phi');
            ccpe_E     = model.getProp(state, 'ccpe_E');
            
            elyte = model.getSubModel('elyte');
            ne    = model.getSubModel('ne');
            pe    = model.getSubModel('pe');
            ccne  = model.getSubModel('ccne');
            ccpe  = model.getSubModel('ccpe');
            
            sp_Li_c  = elyte_c_Li ./ elyte.eps;
            sp_PF6_c = sp_Li_c;

            ne.am.Li.cs = ne.am.Li.cseps ./ ne.am.eps;
            pe.am.Li.cs = pe.am.Li.cseps ./ pe.am.eps;

            %% Update electrolyte physicochemical and transport properties
            state = elyte.update(state);
            
            elyte.kappaeff = elyte.kappa .* elyte.eps .^1.5;
            elyte.sp.Li.Deff = elyte.sp.Li.D .* elyte.eps .^1.5;

            ne.am.update()
            ne.am.Li.Deff = ne.am.Li.D .* ne.am.eps.^1.5;
            ne.sigmaeff = ne.am.sigma .* ne.am.eps.^1.5;

            pe.am.update()
            pe.am.Li.Deff = pe.am.Li.D .* pe.am.eps.^1.5;
            pe.sigmaeff = pe.am.sigma .* pe.am.eps.^1.5;

            obj.ccne.sigmaeff = obj.ccne.am.sigma .* obj.ccne.am.eps.^1.5;
            obj.ccpe.sigmaeff = obj.ccpe.am.sigma .* obj.ccpe.am.eps.^1.5;

            %% Ionic current density
            % Ionic current density in the liquid
            %   Ionic current density due to the chemical potential gradient
            N = elyte.N;
            ncomp = elyte.ncomp;

            jchems = cell(ncomp, 1);
            for i = 1 : ncomp
                jchems{i} = zeros(N + 1, 1);
                coeff =  elyte.kappaeff .* elyte.ion.tvec{i} .* elyte.ion.dmudc{i} ./ ...
                                     (elyte.ion.zvec{i}.*obj.con.F);
                jchems{i} = elyte.operators.harmFace(coeff).* elyte.operators.Grad(elyte.ion.cvec{i});

            end
            elyte.jchem = jchems{1};
            for i = 2 : ncomp
                elyte.jchem = elyte.jchem + jchems{i};
            end
            %   Ionic current density due to the electrochemical potential gradient
            elyte.j = elyte.operators.harmFace(elyte.kappaeff).*(-1).*elyte.operators.Grad(elyte.phi) - elyte.jchem;

            %% Electric current density
            % Active material NE

            ccne = obj.ccne;
            ne = ne;

            ne.j = ne.operators.harmFace(ne.sigmaeff) .* (-1) .* ne.operators.Grad(ne.am.phi);
            obj.ccne.j = ccne.operators.harmFace(ccne.sigmaeff) .* (-1) .* ccne.operators.Grad(ccne.am.phi);

            % Add current transfers between ccne collector and ne material. They correspond to flux continuity
            coupterm = obj.getCoupTerm('ccne-ne');
            face_ccne = coupterm.couplingfaces(:, 1);
            face_ne = coupterm.couplingfaces(:, 2);
            [tne, bccell_ne] = ne.operators.harmFaceBC(ne.sigmaeff, face_ne);
            [tccne, bccell_ccne] = ccne.operators.harmFaceBC(ccne.sigmaeff, face_ccne);
            bcphi_ne = ne.am.phi(bccell_ne);
            bcphi_ccne = ccne.am.phi(bccell_ccne);

            ne.j_bcsource = ne.am.phi*0.0; %NB hack to initialize zero ad
            obj.ccne.j_bcsource = ccne.am.phi*0.0; %NB hack to initialize zero ad

            trans = 1./(1./tne + 1./tccne);
            crosscurrent = trans.*(bcphi_ccne - bcphi_ne);
            ne.j_bcsource(bccell_ne) = crosscurrent;
            obj.ccne.j_bcsource(bccell_ccne) = -crosscurrent;

            % We impose the boundary condition at chosen boundary cells of the ne current collector
            coupterm = obj.getCoupTerm('bc-ccne');
            faces = coupterm.couplingfaces;
            bcval = zeros(numel(faces), 1);
            [tccne, cells] = obj.ccne.operators.harmFaceBC(obj.ccne.sigmaeff, faces);
            obj.ccne.j_bcsource(cells) = obj.ccne.j_bcsource(cells) + tccne.*(bcval - obj.ccne.am.phi(cells));

            % Active material PE and current collector

            ccpe = obj.ccpe;
            pe = pe;

            pe.j =  pe.operators.harmFace(pe.sigmaeff) .* (-1) .* pe.operators.Grad(pe.am.phi);
            obj.ccpe.j =  ccpe.operators.harmFace(ccpe.sigmaeff) .* (-1) .* ccpe.operators.Grad(ccpe.am.phi);

            % Add current transfers between ccpe collector and pe material. They correspond to flux continuity
            coupterm = obj.getCoupTerm('ccpe-pe');
            face_ccpe = coupterm.couplingfaces(:, 1);
            face_pe = coupterm.couplingfaces(:, 2);
            [tpe, bccell_pe] = pe.operators.harmFaceBC(pe.sigmaeff, face_pe);
            [tccpe, bccell_ccpe] = ccpe.operators.harmFaceBC(ccpe.sigmaeff, face_ccpe);
            bcphi_pe = pe.am.phi(bccell_pe);
            bcphi_ccpe = ccpe.am.phi(bccell_ccpe);

            pe.j_bcsource   = pe.am.phi*0.0; %NB hack to initialize zero ad
            obj.ccpe.j_bcsource = ccpe.am.phi*0.0; %NB hack to initialize zero ad

            trans = 1./(1./tpe + 1./tccpe);
            crosscurrent = trans.*(bcphi_ccpe - bcphi_pe);
            pe.j_bcsource(bccell_pe) = crosscurrent;
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
            elyte.sp.Li.source = zeros(elyte.N, 1);
            ne.am.Li.source    = zeros(ne.N, 1);
            pe.am.Li.source    = zeros(pe.N, 1);
            ne.am.e.source     = zeros(ne.N, 1);
            pe.am.e.source     = zeros(pe.N, 1);

            if useAD
                elyte.sp.Li.source = adbackend.convertToAD(elyte.sp.Li.source, adsample);
                ne.am.Li.source    = adbackend.convertToAD(ne.am.Li.source, adsample);
                pe.am.Li.source    = adbackend.convertToAD(pe.am.Li.source, adsample);
                ne.am.e.source     = adbackend.convertToAD(ne.am.e.source, adsample);
                pe.am.e.source     = adbackend.convertToAD(pe.am.e.source, adsample);
            end

            %%%%% Set up chemical source terms %%%%%%%%%%%%%%%%%k
            
            %%%%% NE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            coupterm = obj.getCoupTerm('ne-elyte');
            necells = coupterm.couplingcells(:, 1);
            elytecells = coupterm.couplingcells(:, 2);

            % calculate rection rate
            ne.reactBV(elyte.phi(elytecells));

            % Electrolyte NE Li+ source
            elyte.sp.Li.source(elytecells) = +1 .* ne.R;
            
            % Active Material NE Li0 source
            ne.am.Li.source(necells) = -1 .* ne.R;
            
            % Active Material NE current source
            ne.am.e.source(necells) = +1 .* ne.R;

            %%%%% PE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Electrolyte PE Li+ source
            coupterm = obj.getCoupTerm('pe-elyte');
            pecells = coupterm.couplingcells(:, 1);
            elytecells = coupterm.couplingcells(:, 2);
            
            % calculate rection rate
            pe.reactBV(elyte.phi(elytecells));

            % Electrolyte PE Li+ source
            elyte.sp.Li.source(elytecells) = -1 .* pe.R;

            % Active Material PE Li0 source
            pe.am.Li.source(pecells) = +1 .* pe.R;

            % Active Material PE current source
            pe.am.e.source(pecells) = -1 .* pe.R;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Diffusion Flux                                           %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Divergence of diffusion mass flux
            % Electrolyte Li+ Diffusion

            x = elyte.sp.Li.c;
            flux = - elyte.sp.Li.Trans.*elyte.operators.Grad(x);
            elyte.sp.Li.divDiff = elyte.operators.Div(flux)./elyte.Grid.cells.volumes;

            x= ne.am.Li.cs;
            flux = - ne.am.Li.Trans.*ne.operators.Grad(x);
            ne.am.Li.divDiff =  ne.operators.Div(flux)./ne.Grid.cells.volumes;

            x = pe.am.Li.cs;
            flux = - pe.am.Li.Trans.*pe.operators.Grad(x);
            pe.am.Li.divDiff = pe.operators.Div(flux)./pe.Grid.cells.volumes;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Migration Flux                                           %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Divergence of the migration mass flux
            %   Electrolyte Li+ Migration
            flux = elyte.sp.Li.t ./ (elyte.sp.Li.z .* obj.con.F) .* elyte.j;

            elyte.sp.Li.divMig = elyte.operators.Div(flux)./elyte.Grid.cells.volumes;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% System of Equations                                      %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %% Liquid electrolyte dissolved ionic species mass continuity %
            %   Electrolyte Li+ Mass Continuity
            elyte.sp.Li.massCont = ( - elyte.sp.Li.divDiff - elyte.sp.Li.divMig + elyte.sp.Li.source ...
                                         - elyte.sp.Li.cepsdot);

            %% Liquid electrolyte charge continuity %%%%%%%%%%%%%%%%%%%%%%%

            elyte.chargeCont = -(elyte.operators.Div( elyte.j)./elyte.Grid.cells.volumes) ./ obj.con.F+ ...
                elyte.sp.Li.source .* elyte.sp.Li.z;

            %% Active material mass continuity %%%%%%%%%%%%%%%%%%%%%%%%%%%
            ne.am.Li.massCont = (-ne.am.Li.divDiff + ne.am.Li.source - ne.am.Li.csepsdot);
            pe.am.Li.massCont = (-pe.am.Li.divDiff + pe.am.Li.source - pe.am.Li.csepsdot);

            %% Active material charge continuity %%%%%%%%%%%%%%%%%%%%%%%%%%%

            ne.am.e.chargeCont = (ne.operators.Div(ne.j) - ne.j_bcsource)./ ...
                ne.Grid.cells.volumes./obj.con.F - ne.am.e.source;
            pe.am.e.chargeCont = (pe.operators.Div(pe.j) - pe.j_bcsource)./ ...
                pe.Grid.cells.volumes./obj.con.F - pe.am.e.source;

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
            obj.soe = vertcat(elyte.sp.Li.massCont, ...
                              elyte.chargeCont    , ...
                              ne.am.Li.massCont   , ...
                              ne.am.e.chargeCont  , ...
                              pe.am.Li.massCont   , ...
                              pe.am.e.chargeCont  , ...
                              obj.ccne.am.e.chargeCont, ...
                              obj.ccpe.am.e.chargeCont, ...
                              control);

        end

    end
    
    
end
