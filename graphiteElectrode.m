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
            model = model.initiateCompositeModel();
            
            model.bin  = ptfe();
            model.sei  = seiAM();
            model.eps  = graphitemodel.eps + model.bin.eps + model.sei.eps;
            model.t    = 10e-6;
            
        end

        
        function state = initializeState(model, state)

            state = model.validateState(state);

            graphitemodel = model.getSubModel('graphite');
            state = graphitemodel.initializeState(state);
            OCP   = graphitemodel.getProp(state, 'OCP');
        end
        
        function [namespaces, names] = getModelVarNames(model)
            
            [namespaces1, names1] = getModelVarNames@CompositeModel(model);
            
            names2 = {'eta',  ... % Overpotential,        [V]
                      'j'  ,  ... % Current density,      [A/m2]
                      'R'  ...    % Reaction Rate,
                     };
            namespaces2 = model.assignCurrentNameSpace(names2);
            
            namespaces = horzcat(namespaces1, namespaces2);
            names = horzcat(names1, names2);
            
        end

        function [namespaces, names] = getVarNames(model)
            [namespaces1, names1] = model.getVarNames@CompositeModel();
            namespaces2 = repmat({{}}, 1, 3);
            names2 = {'T', 'SOC', 'phielyte'};

            namespaces = horzcat(namespaces1, namespaces2);
            names = horzcat(names1, names2);
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
        
    end
end

                        
