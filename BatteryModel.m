classdef BatteryModel < CompositeModel

    properties
        fv
    end
    
    methods
        
        function model = BatteryModel(varargin)
            
            model = model@CompositeModel('battery');
            
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
            model.G = G;
            
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
            ccpe.names = {ccpe.names{:}, 'E'};
            model = model.setSubModel(ccpe, 'ccpe');
            
            model.names = {'T', 'SOC'};
            model.aliases = {{'phielyte', VarName({'elyte'}, 'c_Li')}};
            
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
            pvarnames = model.getModelPrimaryVarNames();
            for ind = 1 : numel(pvarnames);
                y0 = [y0; model.getProp(state, pvarnames{ind})];
            end
            
            yp0  = zeros(length(y0), 1);
            
        end

        
        function soe = dynamicBuildSOE(model, t, y, yp, varargin)

            opt = struct('useAD', false);
            opt = merge_options(opt, varargin{:});
            useAD = opt.useAD;

            fv = model.fv;

            if useAD
                adbackend = model.AutoDiffBackend();
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
            model.elyte.sp.Li.cepsdot = yp(fv.getSlot('elyte_Li'));
            model.ne.am.Li.csepsdot   = yp(fv.getSlot('ne_Li'));
            model.pe.am.Li.csepsdot   = yp(fv.getSlot('pe_Li'));

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

            model.ccne.sigmaeff = model.ccne.am.sigma .* model.ccne.am.eps.^1.5;
            model.ccpe.sigmaeff = model.ccpe.am.sigma .* model.ccpe.am.eps.^1.5;

            state = elyte.updateChemicalFluxes(state);
            elyte_j = model.getProp(state, 'j');

            %% Electric current density
            % Active material NE

            state = ccne.updateFlux(model, state);
            state = ne.updateFlux(model, state);

            % Add current transfers between ccne collector and ne material. They correspond to flux continuity
            ne_phi = model.getProp(state, {'ne', 'phi'});
            ccne_phi = model.getProp(state, {'ccne', 'phi'});

            coupterm = model.getCoupTerm('ccne-ne');
            face_ccne = coupterm.couplingfaces(:, 1);
            face_ne = coupterm.couplingfaces(:, 2);
            [tne, bccell_ne] = ne.operators.harmFaceBC(ne.sigmaeff, face_ne);
            [tccne, bccell_ccne] = ccne.operators.harmFaceBC(ccne.sigmaeff, face_ccne);
            
            bcphi_ne = ne_phi(bccell_ne);
            bcphi_ccne = ccne_phi(bccell_ccne);

            ne_j_bcsource = ne_phi*0.0; %NB hack to initialize zero ad
            ccne_j_bcsource = ccne_phi*0.0; %NB hack to initialize zero ad

            trans = 1./(1./tne + 1./tccne);
            crosscurrent = trans.*(bcphi_ccne - bcphi_ne);
            ne_j_bcsource(bccell_ne) = crosscurrent;
            ccne_j_bcsource(bccell_ccne) = -crosscurrent;

            % We impose the boundary condition at chosen boundary cells of the ne current collector
            coupterm = model.getCoupTerm('bc-ccne');
            faces = coupterm.couplingfaces;
            bcval = zeros(numel(faces), 1);
            [tccne, cells] = ccne.operators.harmFaceBC(ccne.sigmaeff, faces);
            ccne_j_bcsource(cells) = ccne_j_bcsource(cells) + tccne.*(bcval - ccne_phi(cells));

            % Active material PE and current collector

            state = ccpe.updateFlux(model, state);
            state = pe.updateFlux(model, state);
            

            % Add current transfers between ccpe collector and pe material. They correspond to flux continuity
            pe_phi = model.getProp(state, {'pe', 'phi'});
            ccpe_phi = model.getProp(state, {'ccpe', 'phi'});

            coupterm = model.getCoupTerm('ccpe-pe');
            face_ccpe = coupterm.couplingfaces(:, 1);
            face_pe = coupterm.couplingfaces(:, 2);
            [tpe, bccell_pe] = pe.operators.harmFaceBC(pe.sigmaeff, face_pe);
            [tccpe, bccell_ccpe] = ccpe.operators.harmFaceBC(ccpe.sigmaeff, face_ccpe);
            bcphi_pe = pe.am.phi(bccell_pe);
            bcphi_ccpe = ccpe.am.phi(bccell_ccpe);

            pe_j_bcsource   = pe_phi*0.0; %NB hack to initialize zero ad
            ccpe_j_bcsource = ccpe_phi*0.0; %NB hack to initialize zero ad

            trans = 1./(1./tpe + 1./tccpe);
            crosscurrent = trans.*(bcphi_ccpe - bcphi_pe);
            pe_j_bcsource(bccell_pe) = crosscurrent;
            ccpe_j_bcsource(bccell_ccpe) = -crosscurrent;

            % We impose the boundary condition at chosen boundary cells of the anode current collector
            coupterm = model.getCoupTerm('bc-ccpe');
            faces = coupterm.couplingfaces;
            bcval = model.ccpe.E;
            [tccpe, cells] = model.ccpe.operators.harmFaceBC(model.ccpe.sigmaeff, faces);
            ccpe_j_bcsource(cells) = ccpe_j_bcsource(cells) + tccpe.*(bcval - ccpe_phi(cells));

            %% Cell voltage
            ccne_E = 0;
            ccpe_E = ccpe.getProp(state, 'E');
            
            U = ccpe_E - ccne_E;
            
            if useAD
                adsample = getSampleAD(y, yp);
                adbackend = model.AutoDiffBackend;
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Source terms for continuity equations                    %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Inititalize vectors
            elyte_Li_source = zeros(elyte.N, 1);
            ne_Li_source    = zeros(ne.N, 1);
            pe_Li_source    = zeros(pe.N, 1);
            ne_e_source     = zeros(ne.N, 1);
            pe_e_source     = zeros(pe.N, 1);

            if useAD
                elyte_Li_source = adbackend.convertToAD(elyte_Li_source, adsample);
                ne_Li_source    = adbackend.convertToAD(ne_Li_source, adsample);
                pe_Li_source    = adbackend.convertToAD(pe_Li_source, adsample);
                ne_e_source     = adbackend.convertToAD(ne_e_source, adsample);
                pe_e_source     = adbackend.convertToAD(pe_e_source, adsample);
            end

            %%%%% Set up chemical source terms %%%%%%%%%%%%%%%%%k
            
            %%%%% NE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            coupterm = model.getCoupTerm('ne-elyte');
            necells = coupterm.couplingcells(:, 1);
            elytecells = coupterm.couplingcells(:, 2);

            % calculate rection rate
            
            error('coupling problem : elyte_phi in ne and in elyte do not match in size');
            state = ne.updateReactBV(model, state);
            ne_R = ne.getProp(state, 'R')

            % Electrolyte NE Li+ source
            elyte_Li_source(elytecells) = ne_R;
            
            % Active Material NE Li0 source
            ne_Li_source(necells) = - ne.R;
            
            % Active Material NE current source
            ne_e_source(necells) = + ne.R;

            %%%%% PE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Electrolyte PE Li+ source
            coupterm = model.getCoupTerm('pe-elyte');
            pecells = coupterm.couplingcells(:, 1);
            elytecells = coupterm.couplingcells(:, 2);
            
            % calculate rection rate
            error('coupling problem : elyte_phi in ne and in elyte do not match in size');
            state = pe.updateReactBV(model, state);
            pe_R = pe.getProp(state, 'R')

            % Electrolyte PE Li+ source
            elyte_Li_source(elytecells) = - pe.R;

            % Active Material PE Li0 source
            pe_Li_source(pecells) = + pe.R;

            % Active Material PE current source
            pe_e_source(pecells) = - pe.R;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Diffusion Flux                                           %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Divergence of diffusion mass flux
            % Electrolyte Li+ Diffusion

            x = elyte.sp.Li.c;
            flux = - elyte.sp.Li.Trans.*elyte.operators.Grad(x);
            elyte_Li_divDiff = elyte.operators.Div(flux)./elyte.G.cells.volumes;

            error('problem between ceps and c');
            x = ne.getProp({'graphite', 'c_Li'});
            flux = - ne.am.Li.Trans.*ne.operators.Grad(x);
            ne_Li_divDiff =  ne.operators.Div(flux)./ne.G.cells.volumes;

            error('problem between ceps and c');
            x = pe.getProp({'nmc111', 'c_Li'});
            flux = - pe.am.Li.Trans.*pe.operators.Grad(x);
            pe_Li_divDiff = pe.operators.Div(flux)./pe.G.cells.volumes;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Migration Flux                                           %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Divergence of the migration mass flux
            %   Electrolyte Li+ Migration
            flux = elyte.sp.Li.t ./ (elyte.sp.Li.z .* model.con.F) .* elyte_j;

            elyte.sp.Li.divMig = elyte.operators.Div(flux)./elyte.G.cells.volumes;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% System of Equations                                      %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %% Liquid electrolyte dissolved ionic species mass continuity %
            %   Electrolyte Li+ Mass Continuity
            elyte.sp.Li.massCont = ( - elyte.sp.Li.divDiff - elyte.sp.Li.divMig + elyte.sp.Li.source ...
                                         - elyte.sp.Li.cepsdot);

            %% Liquid electrolyte charge continuity %%%%%%%%%%%%%%%%%%%%%%%

            elyte.chargeCont = -(elyte.operators.Div( elyte.j)./elyte.G.cells.volumes) ./ model.con.F+ ...
                elyte.sp.Li.source .* elyte.sp.Li.z;

            %% Active material mass continuity %%%%%%%%%%%%%%%%%%%%%%%%%%%
            ne.am.Li.massCont = (-ne.am.Li.divDiff + ne.am.Li.source - ne.am.Li.csepsdot);
            pe.am.Li.massCont = (-pe.am.Li.divDiff + pe.am.Li.source - pe.am.Li.csepsdot);

            %% Active material charge continuity %%%%%%%%%%%%%%%%%%%%%%%%%%%

            ne.am.e.chargeCont = (ne.operators.Div(ne.j) - ne.j_bcsource)./ ...
                ne.G.cells.volumes./model.con.F - ne.am.e.source;
            pe.am.e.chargeCont = (pe.operators.Div(pe.j) - pe.j_bcsource)./ ...
                pe.G.cells.volumes./model.con.F - pe.am.e.source;

            model.ccne.am.e.chargeCont = (model.ccne.operators.Div(model.ccne.j) - model.ccne.j_bcsource)./ ...
                model.ccne.G.cells.volumes./model.con.F;
            model.ccpe.am.e.chargeCont = (model.ccpe.operators.Div(model.ccpe.j) - model.ccpe.j_bcsource)./ ...
                model.ccpe.G.cells.volumes./model.con.F;

            %% control equation
            src = currentSource(t, fv.tUp, fv.tf, model.J);
            coupterm = model.getCoupTerm('bc-ccpe');
            faces = coupterm.couplingfaces;
            bcval = model.ccpe.E;
            [tccpe, cells] = model.ccpe.operators.harmFaceBC(model.ccpe.sigmaeff, faces);
            control = src - sum(tccpe.*(bcval - model.ccpe.am.phi(cells)));

            %% State vector %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            soe = vertcat(elyte.sp.Li.massCont, ...
                          elyte.chargeCont    , ...
                          ne.am.Li.massCont   , ...
                          ne.am.e.chargeCont  , ...
                          pe.am.Li.massCont   , ...
                          pe.am.e.chargeCont  , ...
                          model.ccne.am.e.chargeCont, ...
                          model.ccpe.am.e.chargeCont, ...
                          control);

        end

    end
    
    
end
