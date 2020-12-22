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
        
        % Material properties
        E   % Electric potential,   [V]
        
        bin     % Binder object
        ca      % Conducting additive object
        cei     % Cathode-electrolyte interphase (CEI) object
        elyte   % Liquid electrolyte data structure
        
        % Effective conductivity
        sigmaeff

        
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
            nmc111model = nmc111AM();
            model.SubModels{1} = nmc111model;
            model.SubModelNames{1} = 'nmc111';
            
            model.bin = ptfe();
            model.eps = nmc111model.eps + model.bin.eps;
            model.t = 10e-6;
            
            model.aliases = {{'T', VarName({}, 'T')}, ...
                             {'SOC', VarName({}, 'SOC')}, ...
                             {'phielyte', VarName({}, 'phielyte')}, ...
                            };
            
            model = model.initiateCompositeModel();
            
        end
        
        function state = initializeState(model, state)

            state = model.validateState(state);

            nmc111model = model.getSubModel('nmc111');
            state = nmc111model.initializeState(state);
            OCP   = nmc111model.getProp(state, 'OCP');
        end

        
        function state = reactBV(model, state)
            
            nmc111model = model.getSubModel('nmc111');
            
            state = nmc111model.updateEquilibrium(state);

            T   = model.getProp(state, 'T');
            eta = model.getProp(state, 'eta');
            phiElite = model.getProp(state, 'phielyte');
            
            phi = nmc111model.getProp(state, 'phi');
            OCP = nmc111model.getProp(state, 'OCP');
            k   = nmc111model.getProp(state, 'k');
            
            eta = - (phi - phiElyte - OCP);
                                    
            R = nmc111model.Asp.*butlerVolmer(k.*model.con.F, 0.5, 1, eta, T) ./ (1 .* model.con.F);
            
            state = model.setProp(state, 'R');
            
        end
    end
end

