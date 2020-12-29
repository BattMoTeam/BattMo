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

            warning('setup eps');
            
            %% setup elyte
            nx = sum(nxs); 
            
            submodels = {};
            
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

            model.SubModels = submodels;
            model.hasparent = false;
            
            model = model.initiateCompositeModel();

            ccpe = model.getAssocModel('ccpe');
            ccpe.pnames = {ccpe.pnames{:}, 'E'};
            ccpe.names = {ccpe.names{:}, 'E'};
            model = model.setSubModel('ccpe', ccpe);
            
            model.names = {'T', 'SOC'};

            % setup ne
            
            ne = model.getAssocModel('ne');
            
            fnupdate = @(model, state) model.dispatchValues(state);
            fnmodel = {'..'};
            ne = ne.setPropFunction(PropFunction('T', fnupdate, fnmodel));
            ne = ne.setPropFunction(PropFunction('SOC', fnupdate, fnmodel));
            
            fnupdate = @(model, state) model.updatePhiElyte(state);
            fnmodel = {'..'};
            ne = ne.setPropFunction(PropFunction('phielyte', fnupdate, fnmodel));
            
            model = model.setSubModel('ne', ne);
            
            % setup pe
            
            pe = model.getAssocModel('pe');
            
            fnupdate = @(model, state) model.dispatchValues(state);
            fnmodel = {'..'};
            pe = pe.setPropFunction(PropFunction('T', fnupdate, fnmodel));
            pe = pe.setPropFunction(PropFunction('SOC', fnupdate, fnmodel));
            
            fnupdate = @(model, state) model.updatePhiElyte(state);
            fnmodel = {'..'};
            pe = pe.setPropFunction(PropFunction('phielyte', fnupdate, fnmodel));
            
            model = model.setSubModel('pe', pe);
            
            % setup ccne
            ccne = model.getAssocModel('ccne');
            fnupdate = @(model, state) model.dispatchValues(state);
            fnmodel = {'..'};
            ccne = ccne.setPropFunction(PropFunction('T', fnupdate, fnmodel));
            model = model.setSubModel('ccne', ccne);
            
            % setup ccpe
            ccpe = model.getAssocModel('ccpe');
            fnupdate = @(model, state) model.dispatchValues(state);
            fnmodel = {'..'};
            ccpe = ccpe.setPropFunction(PropFunction('T', fnupdate, fnmodel));
            model = model.setSubModel('ccpe', ccpe);
            
            % setup elyte
            elyte = model.getAssocModel('elyte');
            fnupdate = @(model, state) model.dispatchValues(state);
            fnmodel = {'..'};
            elyte = elyte.setPropFunction(PropFunction('T', fnupdate, fnmodel));
            model = model.setSubModel('elyte', elyte);
            
            model = model.initiateCompositeModel();
            
        end

        function state = initializeState(model, state)
            
            state = model.validateState(state);
            
            % initialize each of the submodels
            for i = 1 : model.nSubModels
                state = model.SubModels{i}.initializeState(state);
            end
            
            ccne = model.getAssocModel('ccne');
            ne = model.getAssocModel('ne');
            
            [OCP_ne, state] = ne.getUpdatedProp(state, {'am', 'OCP'});
            nc = ccne.G.cells.num;
            OCP_ccne = OCP_ne(1)*ones(nc, 1);
            
            state = ccne.setProp(state, 'OCP', OCP_ccne);
            state = ccne.setProp(state, 'phi', OCP_ccne);
            
            ne = model.getAssocModel('ccpe');
            pe = model.getAssocModel('pe');
            ccpe = model.getAssocModel('ccpe');

            [OCP_pe, state] = pe.getUpdatedProp(state, {'am', 'OCP'});
            nc = ccpe.G.cells.num;
            OCP_ccpe = OCP_pe(1)*ones(nc, 1);
            
            state = ccpe.setProp(state, 'OCP', OCP_ccpe);
            state = ccpe.setProp(state, 'phi', OCP_ccpe);
            
            state = ccpe.setProp(state, 'E', OCP_ccpe(1));
            
        end

        
        function state = dispatchValues(model, state)
            
            [T, state] = model.getUpdatedProp(state, 'T');
            [SOC, state] = model.getUpdatedProp(state, 'SOC');
            
            elyte = model.getAssocModel('elyte');
            G = elyte.G;
            Telyte = T(G.mappings.cellmap);
            
            ne = model.getAssocModel('ne');
            G = ne.G;
            Tne = T(G.mappings.cellmap);
            SOCne = SOC(G.mappings.cellmap);
            
            pe = model.getAssocModel('pe');
            G = pe.G;
            Tpe = T(G.mappings.cellmap);
            SOCpe = SOC(G.mappings.cellmap);
            
            ccpe = model.getAssocModel('ccpe');
            G = ccpe.G;
            Tccpe = T(G.mappings.cellmap);
            
            ccne = model.getAssocModel('ccne');
            G = ccne.G;
            Tccne = T(G.mappings.cellmap);
            
            state = model.setProp(state, {'elyte', 'T'}, Telyte);
            state = model.setProp(state, {'ne', 'T'}, Tne);
            state = model.setProp(state, {'pe', 'T'}, Tpe);
            state = model.setProp(state, {'ccpe', 'T'}, Tccpe);
            state = model.setProp(state, {'ccne', 'T'}, Tccne);
            state = model.setProp(state, {'ne', 'SOC'}, SOCne);
            state = model.setProp(state, {'pe', 'SOC'}, SOCpe);
            
        end

        function state = updatePhiElyte(model, state)
            
            [phielyte, state] = model.getUpdatedProp(state, {'elyte', 'phi'});
            
            ne = model.getAssocModel('ne');
            G = ne.G;
            phine = phielyte(G.cellmap);
            
            pe = model.getAssocModel('pe');
            G = pe.G;
            phipe = phielyte(G.cellmap);
            
            state = model.setProp(state, {'ne', 'phielyte'}, phine);
            state = model.setProp(state, {'pe', 'phielyte'}, phipe);
            
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
            model.elyte.Li.cepsdot = yp(fv.getSlot('elyte_Li'));
            model.ne.am.Li.csepsdot   = yp(fv.getSlot('ne_Li'));
            model.pe.am.Li.csepsdot   = yp(fv.getSlot('pe_Li'));

            elyte = model.getAssocModel('elyte');
            ne    = model.getAssocModel('ne');
            pe    = model.getAssocModel('pe');
            ccne  = model.getAssocModel('ccne');
            ccpe  = model.getAssocModel('ccpe');

            elyte_c_Li = elyte.getUpdatedProp(state, 'c_Li');
            elyte_phi  = elyte.getUpdatedProp(state, 'phi');
            ne_Li      = ne.getUpdatedProp(state, {'am', 'Li'});
            ne_phi     = ne.getUpdatedProp(state, {'am', 'phi'});
            pe_Li      = pe.getUpdatedProp(state, {'am', 'Li'});
            pe_phi     = pe.getUpdatedProp(state, {'am', 'phi'});
            ccne_phi   = ccne.getUpdatedProp(state, 'phi');
            ccpe_phi   = ccpe.getUpdatedProp(state, 'phi');
            ccpe_E     = ccpe.getUpdatedProp(state, 'E');
            
            sp_Li_c  = elyte_c_Li ./ elyte.eps;
            sp_PF6_c = sp_Li_c;

            ne.am.Li.cs = ne.am.Li.cseps ./ ne.am.eps;
            pe.am.Li.cs = pe.am.Li.cseps ./ pe.am.eps;

            %% We compute the effective conductivity and diffusion
            
            [D, state] = elyte.getUpdatedProp(state, 'D');
            [sigma, state] = elyte.getUpdatedProp(state, 'D');
            elyte_Deff = D .* elyte.eps .^1.5;
            elyte_sigmaeff = sigma .* elyte.eps .^1.5;            

            [D, state] = ne.getUpdatedProp(state, 'D');
            [sigma, state] = ne.getUpdatedProp(state, 'D');
            ne_Deff = D .* ne.eps .^1.5;
            ne_sigmaeff = sigma .* ne.eps .^1.5;            

            [D, state] = pe.getUpdatedProp(state, 'D');
            [sigma, state] = pe.getUpdatedProp(state, 'D');
            pe_Deff = D .* pe.eps .^1.5;
            pe_sigmaeff = sigma .* pe.eps .^1.5;            

            ccne_sigmaeff = ccne.sigmaeff;
            ccpe_sigmaeff = ccpe.sigmaeff;
            
            %% We setup the current transfers between ccne collector and ne material with ne_j_bcsource and ccne_j_bcsource
            
            ne_j_bcsource = ne_phi*0.0; %NB hack to initialize zero ad
            ccne_j_bcsource = ccne_phi*0.0; %NB hack to initialize zero ad

            coupterm = model.getCoupTerm('ccne-ne');
            face_ccne = coupterm.couplingfaces(:, 1);
            face_ne = coupterm.couplingfaces(:, 2);
            [tne, bccell_ne] = ne.operators.harmFaceBC(ne.sigmaeff, face_ne);
            [tccne, bccell_ccne] = ccne.operators.harmFaceBC(ccne_sigmaeff, face_ccne);
            
            bcphi_ne = ne_phi(bccell_ne);
            bcphi_ccne = ccne_phi(bccell_ccne);

            trans = 1./(1./tne + 1./tccne);
            crosscurrent = trans.*(bcphi_ccne - bcphi_ne);
            ne_j_bcsource(bccell_ne) = crosscurrent;
            ccne_j_bcsource(bccell_ccne) = -crosscurrent;

            %% We impose the boundary condition at chosen boundary cells of the ne current collector by updating ccne_j_bcsource
            
            coupterm = model.getCoupTerm('bc-ccne');
            faces = coupterm.couplingfaces;
            bcval = zeros(numel(faces), 1);
            [tccne, cells] = ccne.operators.harmFaceBC(ccne.sigmaeff, faces);
            ccne_j_bcsource(cells) = ccne_j_bcsource(cells) + tccne.*(bcval - ccne_phi(cells));

            %% We setup the current transfers between ccpe collector and pe material with pe_j_bcsource and ccpe_j_bcsource

            pe_j_bcsource   = pe_phi*0.0; %NB hack to initialize zero ad
            ccpe_j_bcsource = ccpe_phi*0.0; %NB hack to initialize zero ad

            coupterm = model.getCoupTerm('ccpe-pe');
            face_ccpe = coupterm.couplingfaces(:, 1);
            face_pe = coupterm.couplingfaces(:, 2);
            [tpe, bccell_pe] = pe.operators.harmFaceBC(pe.sigmaeff, face_pe);
            [tccpe, bccell_ccpe] = ccpe.operators.harmFaceBC(ccpe.sigmaeff, face_ccpe);
            bcphi_pe = pe.am.phi(bccell_pe);
            bcphi_ccpe = ccpe.am.phi(bccell_ccpe);

            trans = 1./(1./tpe + 1./tccpe);
            crosscurrent = trans.*(bcphi_ccpe - bcphi_pe);
            pe_j_bcsource(bccell_pe) = crosscurrent;
            ccpe_j_bcsource(bccell_ccpe) = -crosscurrent;

            %% We impose the boundary condition at chosen boundary cells of the anode current collector by updating ccpe_j_bcsource
            
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

            % We setup the exchange terms between the electrolyte and the electrodes for Li due to chemical reacions:
            % elyte_Li_source, ne_Li_source, pe_Li_source.
            %
            % We setup also the electron production in the electrod: ne_e_source, pe_e_source.
            %
            
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
            [ne_R, state] = ne.getUpdatedProp(state, 'R')

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
            [pe_R, state] = pe.getUpdatedProp(state, 'R')

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
            
            trans = elyte.operators.harmFace(elyte_Deff);
            flux = - trans.*elyte.operators.Grad(x);
            elyte_Li_divDiff = elyte.operators.Div(flux)./elyte.G.cells.volumes;

            warning('problem between ceps and c');
            x = ne.getUpdatedProp(state, 'c_Li');
            trans = ne.operators.harmFace(ne_Deff);
            flux = - trans.*ne.operators.Grad(x);
            ne_Li_divDiff = ne.operators.Div(flux)./ne.G.cells.volumes;
            
            warning('problem between ceps and c');
            x = pe.getUpdatedProp(state, 'c_Li');
            trans = pe.operators.harmFace(pe_Deff);
            flux = - trans.*pe.operators.Grad(x);
            pe_Li_divDiff = pe.operators.Div(flux)./pe.G.cells.volumes;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Migration Flux                                           %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Divergence of the migration mass flux for the electrolyte Li+ Migration
            
            [elyte_j, state] = elyte.getUpdatedProp(state, 'j');
            flux = elyte.sp.t ./ (elyte.sp.z .* model.con.F) .* elyte_j;
            elyte_Li_divMig = elyte.operators.Div(flux)./elyte.G.cells.volumes;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% System of Equations                                      %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %% Liquid electrolyte dissolved ionic species mass continuity %
            %   Electrolyte Li+ Mass Continuity
            elyte_Li_massCont = (-elyte_Li_divDiff - elyte_Li_divMig + elyte_Li_source - elyte_Li_cepsdot);

            %% Liquid electrolyte charge continuity %%%%%%%%%%%%%%%%%%%%%%%

            elyte_chargeCont = -(elyte.operators.Div(elyte_j)./elyte.G.cells.volumes)./model.con.F + ...
                elyte_Li_source.*elyte_sp_z;

            %% Electrode Active material mass continuity %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            ne_Li_massCont = (-ne_Li_divDiff + ne_Li_source - ne_Li_csepsdot);
            pe_Li_massCont = (-pe_Li_divDiff + pe_Li_source - pe_Li_csepsdot);

            %% Electode Active material charge continuity %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            [ne_j, state] = ne.getUpdatedProp(state, 'j');
            [pe_j, state] = pe.getUpdatedProp(state, 'j');
            
            ne_e_chargeCont = (ne.operators.Div(ne_j) - ne_j_bcsource)./ ne.G.cells.volumes./model.con.F - ...
                ne_e_source;
            pe_e_chargeCont = (pe.operators.Div(pe_j) - pe_j_bcsource)./ pe.G.cells.volumes./model.con.F - ...
                pe_e_source;

            %% Collector charge continuity
            
            [ccne_j, state] = ccne.getUpdatedProp(state, 'j');
            [ccpe_j, state] = ccpe.getUpdatedProp(state, 'j');
            
            ccne_e_chargeCont = (ccne.operators.Div(ccne_j) - ccne_j_bcsource)./ ccne.G.cells.volumes./model.con.F;
            ccpe_e_chargeCont = (ccpe.operators.Div(ccpe_j) - ccpe.j_bcsource)./ ccpe.G.cells.volumes./model.con.F;_
            
            %% control equation
            
            src = currentSource(t, fv.tUp, fv.tf, model.J);
            coupterm = model.getCoupTerm('bc-ccpe');
            faces = coupterm.couplingfaces;
            bcval = model.ccpe.E;
            [tccpe, cells] = model.ccpe.operators.harmFaceBC(ccpe_sigmaeff, faces);
            control = src - sum(tccpe.*(bcval - ccpe_phi(cells)));

            %% Governing equations 
            
            soe = vertcat(elyte_Li_massCont, ...
                          elyte_chargeCont , ...
                          ne_Li_massCont   , ...
                          ne_e_chargeCont  , ...
                          pe_Li_massCont   , ...
                          pe_e_chargeCont  , ...
                          ccne_e_chargeCont, ...
                          ccpe_e_chargeCont, ...
                          control);

        end

    end
    
    
end
