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
                
    end
    
    methods
        
        function model = graphiteElectrode(name, G, cells)
            
            model = model@CompositeModel(name);
            model.G = genSubGrid(G, cells);
            
            % setup graphite submodel
            graphitemodel = graphiteAM();
            model.SubModels{1} = graphitemodel;
            model.SubModelNames{1} = 'graphite';
            
            model.bin  = ptfe();
            model.sei  = seiAM();
            model.eps  = graphitemodel.eps + model.bin.eps + model.sei.eps;
            model.t    = 10e-6;
            
            % setup sigmaeff
            eps = graphitemodel.eps;
            sigma = graphitemodel.sigma;
            model.sigmaeff = sigma .* eps.^1.5;            
            
            % state variables
            names = {'eta',  ... % Overpotential,        [V]
                     'j'  ,  ... % Current density,      [A/m2]
                     'R'  ...    % Reaction Rate,
                    };
            model.names = names;
            
            model.aliases = {{'T', VarName({}, 'T')}, ...
                             {'SOC', VarName({}, 'SOC')}, ...
                             {'phielyte', VarName({}, 'phielyte')}, ...
                            };
            
            model = model.initiateCompositeModel();
        end

        
        function state = initializeState(model, state)

            state = model.validateState(state);

            graphitemodel = model.getSubModel('graphite');
            state = graphitemodel.initializeState(state);
            OCP   = graphitemodel.getProp(state, 'OCP');
        end
        
        function state = reactBV(model, state)
            
            graphitemodel = model.getSubModel('graphite');
            
            state = graphitemodel.updateEquilibrium(state);

            T   = model.getProp(state, 'T');
            eta = model.getProp(state, 'eta');
            phiElite = model.getProp(state, 'phielyte');
            
            phi = graphitemodel.getProp(state, 'phi');
            OCP = graphitemodel.getProp(state, 'OCP');
            k   = graphitemodel.getProp(state, 'k');
            
            eta = (phi - phiElyte - OCP);
                                    
            R = graphitemodel.Asp.*butlerVolmer(k.*model.con.F, 0.5, 1, eta, T) ./ (1 .* model.con.F);
            
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

                        
