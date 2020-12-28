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
            
            % state variables
            names = {'eta',  ... % Overpotential,        [V]
                     'j',  ... % Current density,      [A/m2]
                     'R'... % Reaction Rate,
                    };
        
            model.names = names;
            
            % setup nmc111 submodel
            amModel = nmc111AM();
            model.SubModels{1} = amModel;
            model.SubModelNames{1} = amModel.getModelName;
            
            model.bin = ptfe();
            model.eps = amModel.eps + model.bin.eps;
            model.t = 10e-6;

            % setup sigmaeff
            eps = amModel.eps;
            sigma = amModel.sigma;
            model.sigmaeff = sigma .* eps.^1.5;
                        
            model.aliases = {{'T', VarName({}, 'T')}, ...
                             {'SOC', VarName({}, 'SOC')}, ...
                             {'phielyte', VarName({}, 'phielyte')}, ...
                            };
            
            model = model.initiateCompositeModel();
            
        end
        
        function state = initializeState(model, state)

            state = model.validateState(state);

            amModel = model.getAssocModel('nmc111');
            state = amModel.initializeState(state);
            OCP   = amModel.getProp(state, 'OCP');
        end

        
        function state = updateReactBV(model, state)
            
            amModel = model.getAssocModel('nmc111');
            
            state = amModel.updateEquilibrium(state);

            T   = model.getProp(state, 'T');
            eta = model.getProp(state, 'eta');
            phiElite = model.getProp(state, 'phielyte');
            
            phi = amModel.getProp(state, 'phi');
            OCP = amModel.getProp(state, 'OCP');
            k   = amModel.getProp(state, 'k');
            
            eta = - (phi - phiElyte - OCP);
                                    
            R = amModel.Asp.*butlerVolmer(k.*model.con.F, 0.5, 1, eta, T) ./ (1 .* model.con.F);
            
            state = model.setProp(state, 'R');
            
        end
        
        function state = updateFlux(model, state)
            
            phi = model.getProp(state, 'phi');
            op = model.operators;
            sigmaeff = model.sigmaeff;

            j = - op.harmFace(sigmaeff).*op.Grad(phi); 
            
            state = model.setProp(state, 'j', j);
            
        end
        
    end
end

