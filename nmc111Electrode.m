classdef nmc111Electrode < CompositeModel
    % NMC111ELECTRODE Summary of this class goes here
    % Detailed explanation goes here
    
    properties
        
        % Physical constants
        con = physicalConstants();
        
        % Design properties
        t       % Thickness,        [m]
        eps     % Volume fraction,  [-]
        void    % Porosity,         [-]
        sigmaeff % effective sigma
        
        % Material properties
        E   % Electric potential,   [V]
        
        bin     % Binder object
        ca      % Conducting additive object
        cei     % Cathode-electrolyte interphase (CEI) object
        elyte   % Liquid electrolyte data structure
        
    end
    
    methods 
        
        function model = nmc111Electrode(name, G, cells)
            
            model = model@CompositeModel(name);
            model.G = genSubGrid(G, cells);
            nc = model.G.cells.num;
            
            % setup operators
            model.operators = localSetupOperators(model.G);

            % setup nmc111 submodel
            ammodel = nmc111AM('am');
            ammodel = ammodel.setAlias({'T', VarName({'..'}, 'T')});
            ammodel = ammodel.setAlias({'SOC', VarName({'..'}, 'SOC')});

            model.SubModels{1} = ammodel;
            
            model.bin = ptfe();
            model.eps = (ammodel.eps + model.bin.eps)*ones(nc, 1);
            model.void = 1 - model.eps;
            model.t = 10e-6;

            % setup sigmaeff
            eps = (ammodel.eps)*ones(nc, 1);
            sigma = ammodel.sigma;
            model.sigmaeff = sigma .* eps.^1.5;

            % state variables
            names = {'R'         , ... % Reaction Rate   ,
                     'jBcSource' , ...
                     'LiSource'  , ...
                     'eSource'   , ...
                     'chargeCont', ...
                     'LiFlux'    , ...
                     'T'         , ... % temperature
                     'SOC'       , ...
                     'phielyte'};
        
            model.names = names;
            
            %% setup Update property functions
            
            propfunctions = {};

            updatefn = @(model, state) model.updateQuantities(state);
            names = {'chargeCont', 'LiFlux'};
            
            for ind = 1 : numel(names)
                name = names{ind};
                propfunction = PropFunction(name, updatefn, '.');
                propfunctions{end + 1} = propfunction;
            end
            
            name = 'R';
            updatefn = @(model, state) model.updateReactionRate(state);
            propfunction = PropFunction(name, updatefn, '.');
            propfunctions{end + 1} = propfunction;
            
            model.propfunctions = propfunctions;

            model = model.initiateCompositeModel();
            
        end
        
        function state = initializeState(model, state)
           % Used only in debugging for the moment
            state = model.initiateState(state);

            ammodel = model.getAssocModel('am');
            state = ammodel.initializeState(state);
        end

        
        function state = updateReactionRate(model, state)
            

            [T, state]   = model.getUpdatedProp(state, 'T');
            [phiElyte, state] = model.getUpdatedProp(state, 'phielyte');

            ammodel = model.getAssocModel('am');

            [phi, state] = ammodel.getUpdatedProp(state, 'phi');
            [OCP, state] = ammodel.getUpdatedProp(state, 'OCP');
            [k, state]   = ammodel.getUpdatedProp(state, 'k');
            
            eta = -(phi - phiElyte - OCP);
                                    
            R = ammodel.Asp.*butlerVolmer(k.*model.con.F, 0.5, 1, eta, T) ./ (1 .* model.con.F);
            
            state = model.setProp(state, 'R', R);
            
        end
        
        function state = updateQuantities(model, state)
            
            op = model.operators;
            [phi, state] = model.getUpdatedProp(state, {'am', 'phi'});
            op = model.operators;
            
            sigmaeff = model.sigmaeff;
                                
            j = - op.harmFace(sigmaeff).*op.Grad(phi); 

            [jBcSource, state] = model.getUpdatedProp(state, 'jBcSource');
            [eSource, state]   = model.getUpdatedProp(state, 'eSource');
            
            chargeCont = (op.Div(j) - jBcSource)./ model.G.cells.volumes./model.con.F - eSource;

            [D, state]   = model.getUpdatedProp(state, {'am', 'D'});
            [cLi, state] = model.getUpdatedProp(state, {'am', 'Li'});
            
            Deff = D .* model.eps .^1.5;
            
            trans = op.harmFace(Deff);
            flux = - trans.*op.Grad(cLi);
            
            state = model.setProp(state, 'chargeCont', chargeCont);
            state = model.setProp(state, 'LiFlux', flux);
            
            
        end
        


    end
end

