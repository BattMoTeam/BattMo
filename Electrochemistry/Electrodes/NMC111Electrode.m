classdef NMC111Electrode < CompositeModel
    % NMC111ELECTRODE Electrode class with NMC111 active material
    % Detailed explanation goes here

    properties
        
        % Physical constants
        constants = PhysicalConstants();
        
        % Design properties
        thickness       % Thickness,        [m]
        volumeFraction  % Volume fraction,  [-]
        porosity        % Porosity,         [-]
        mass            % Mass,             [kg]
        volume          % Volume,           [m3]
        area            % Area,             [m2]
        
        % Transport properties
        effectiveElectronicConductivity % effective electronicConductivity
        
        % Material properties
        electricPotential   % Electric potential,   [V]
        
        % submodels
        ActiveMaterial      % Active material
        Electrolyte         % Electrolyte
        Binder              % Binder object
        ConductingAdditive  % Conducting additive object
        CEI                 % Cathode-electrolyte interphase (CEI) object
        
    end
    
    methods 
        
        function model = NMC111Electrode(name, G, cells)
            
            model = model@CompositeModel(name);
            model.G = genSubGrid(G, cells);
            nc = model.G.cells.num;
            
            % setup operators
            model.operators = localSetupOperators(model.G);

            % setup nmc111 submodel
            ActiveMaterialmodel = NMC111('ActiveMaterial');
            ActiveMaterialmodel = ActiveMaterialmodel.setAlias({'T', VarName({'..'}, 'T')});
            ActiveMaterialmodel = ActiveMaterialmodel.setAlias({'SOC', VarName({'..'}, 'SOC')});

            model.SubModels{1} = ActiveMaterialmodel;
            
            % model.bin = ptfe();
            % model.volumeFraction = (ActiveMaterialmodel.volumeFraction + model.bin.volumeFraction)*ones(nc, 1);
            model.volumeFraction = ActiveMaterialmodel.volumeFraction*ones(nc, 1);
            model.porosity = 1 - model.volumeFraction;
            model.thickness = 10e-6;

            % setup effectiveElectronicConductivity
            eps = (ActiveMaterialmodel.volumeFraction)*ones(nc, 1);
            electronicConductivity = ActiveMaterialmodel.electronicConductivity;
            model.effectiveElectronicConductivity = electronicConductivity .* eps.^1.5;

            % state variables
            names = {'R'         , ... % Reaction Rate   ,
                     'jBcSource' , ...
                     'LiSource'  , ...
                     'eSource'   , ...
                     'chargeCont', ...
                     'massCont'  , ...
                     'LiFlux'    , ...
                     'T'         , ... % temperature
                     'SOC'       , ...
                     'phiElectrolyte'};
        
            model.names = names;
            model = model.setupVarDims();
            
            %% setup Update property functions
            
            names = {'chargeCont', 'LiFlux'};
            updatefn = @(model, state) model.updateQuantities(state);
            inputnames = {'jBcSource', 'eSource'};
            for ind = 1 : numel(names)
                name = names{ind};
                model = model.addPropFunction(name, updatefn, inputnames, {'.'});
            end
            
            name = 'R';
            updatefn = @(model, state) model.updateReactionRate(state);
            inputnames = {'phiElectrolyte', {'ActiveMaterial', 'T'}, {'ActiveMaterial', 'phi'}, {'ActiveMaterial', 'OCP'}, {'ActiveMaterial', 'k'}};
            model = model.addPropFunction(name, updatefn, inputnames, {'.'});

            model = model.initiateCompositeModel();
            
            model.ActiveMaterial = model.getAssocModel('ActiveMaterial');
        
        end
        
        
        function state = updateReactionRate(model, state)
            

            T = state.ActiveMaterial.T;
            phiElyte = state.phiElectrolyte;
            
            ActiveMaterialmodel = model.ActiveMaterial;
                        
            state.ActiveMaterial = ActiveMaterialmodel.updateQuantities(state.ActiveMaterial);
            phi = state.ActiveMaterial.phi;
            OCP = state.ActiveMaterial.OCP;
            k = state.ActiveMaterial.k;

            eta = -(phi - phiElyte - OCP);
            state.eta = eta;                        
            R = ActiveMaterialmodel.volumetricSurfaceArea.*ButlerVolmerEquation(k.*model.constants.F, 0.5, 1, eta, T) ./ (1 .* model.constants.F);
            
            state.R = R;
            
        end
        
        function state = updateQuantities(model, state)
            
            phi = state.ActiveMaterial.phi;

            op = model.operators;
            
            effectiveElectronicConductivity = model.effectiveElectronicConductivity;
                                
            j = - op.harmFace(effectiveElectronicConductivity).*op.Grad(phi); 
                       
            % We assume the source have been computed before this function is run
            jBcSource = state.jBcSource;
            eSource   = state.eSource;
            
            chargeCont = (op.Div(j) - jBcSource)./ model.G.cells.volumes./model.constants.F - eSource;

            % We assume the following quantities have been updated in ActiveMaterial model
            D =  state.ActiveMaterial.D;
            cLi =  state.ActiveMaterial.Li;
            
            Deff = D .* model.volumeFraction .^1.5;
            
            trans = op.harmFace(Deff);
            flux = - trans.*op.Grad(cLi);
            
            state.chargeCont = chargeCont;
            state.LiFlux = flux;
            
        end

    end
end

