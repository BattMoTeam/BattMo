classdef graphiteElectrode < CompositeModel
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        % Physical constants
        con = physicalConstants();
        
        % Design properties
        t       % Thickness,        [m]
        eps     % Volume fraction,  [-]
        void    % Porosity,         [-]
        
        % Material properties
        bin     % Binder object
        ca      % Conducting additive object
        sei     % Solid-electrolyte interphase (SEI) object
        elyte   % Liquid electrolyte data structure
        E       % Electric potential,   [V]
        
        % Effective conductivity
        sigmaeff

        % submodels
        am
        
    end
    
    methods
        
        function model = graphiteElectrode(name, G, cells)
            
            model = model@CompositeModel(name);
            model.G = genSubGrid(G, cells);
            nc = model.G.cells.num;
            
            % setup operators
            model.operators = localSetupOperators(model.G);
            
            % setup graphite submodel
            ammodel = graphiteAM('am');
            ammodel = ammodel.setAlias({'T', VarName({'..'}, 'T')});
            ammodel = ammodel.setAlias({'SOC', VarName({'..'}, 'SOC')});

            model.SubModels{1} = ammodel;
            
            model.bin  = ptfe();
            model.sei  = seiAM();
            model.eps  = (ammodel.eps + model.bin.eps + model.sei.eps)*ones(nc, 1);
            model.void = 1 - model.eps;
            model.t    = 10e-6;
            
            % setup sigmaeff
            eps = (ammodel.eps)*ones(nc, 1);
            sigma = ammodel.sigma;
            model.sigmaeff = sigma .* eps.^1.5;            
            
            % state variables
            names = {'R'         , ...
                     'jBcSource' , ...
                     'LiSource'  , ...
                     'eSource'   , ...
                     'chargeCont', ...
                     'massCont'  , ...
                     'LiFlux'    , ...
                     'T'         , ...
                     'SOC'       , ...
                     'phielyte'};
            model.names = names;
            model = model.setupVarDims();
            
            %% setup Update property functions
            
            propfunctions = {};

            updatefn = @(model, state) model.updateQuantities(state);
            names = {'chargeCont', 'LiFlux'};
            inputnames = {'jBcSource', 'eSource'};            
            for ind = 1 : numel(names)
                name = names{ind};
                model = model.addPropFunction(name, updatefn, inputnames, {'.'});
            end
            
            name = 'R';
            updatefn = @(model, state) model.updateReactionRate(state);
            inputnames = {'phielyte', {'am', 'T'}, {'am', 'phi'}, {'am', 'OCP'}, {'am', 'k'}};
            model = model.addPropFunction(name, updatefn, inputnames, {'.'});
            
            model = model.initiateCompositeModel();
        
            model.am = model.getAssocModel('am');
            
        end
        
        function state = updateReactionRate(model, state)
             
            T = state.T;
            phiElyte = state.phielyte;
            
            ammodel = model.am;
                        
            state.am = ammodel.updateQuantities(state.am);
            phi = state.am.phi;
            OCP = state.am.OCP;
            k = state.am.k;
            
            eta = (phi - phiElyte - OCP);
            
            R = ammodel.Asp.*butlerVolmer(k.*model.con.F, 0.5, 1, eta, T) ./ (1 .* model.con.F);
             
            state.R = R;
            
        end
        
        
        function state = updateQuantities(model, state)
            
            phi = state.am.phi;

            op = model.operators;
            
            sigmaeff = model.sigmaeff;
                                
            j = - op.harmFace(sigmaeff).*op.Grad(phi); 
                       
            % We assume the source have been computed before this function is run
            jBcSource = state.jBcSource;
            eSource   = state.eSource;
            
            chargeCont = (op.Div(j) - jBcSource)./ model.G.cells.volumes./model.con.F - eSource;

            % We assume the following quantities have been updated in am model
            D =  state.am.D;
            cLi =  state.am.Li;
            
            Deff = D .* model.eps .^1.5;
            
            trans = op.harmFace(Deff);
            flux = - trans.*op.Grad(cLi);
            
            state.chargeCont = chargeCont;
            state.LiFlux = flux;
            
        end        
    end
end

                        
