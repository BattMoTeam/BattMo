classdef Battery2 < PhysicalModel

    properties
        
        con = PhysicalConstants();

        % Temperature and SOC
        % for the moment here, for convenience. Will be moved
        T
        SOC

        % Input current
        J
        % Voltage cut
        Ucut

        % Components
        Electrolyte
        NegativeElectrode
        PositiveElectrode
    
        couplingTerms
        
    end
    
    methods
        
        function model = Battery2(paramobj)
            
            model = model@PhysicalModel([]);
            
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks',true);
            
            %% Setup the model using the input parameters
            fdnames = {'G', ...
                       'couplingTerms', ...
                       'T'  , ...
                       'SOC', ...
                       'J'  , ...
                       'Ucut'};
            
            model = dispatchParams(model, paramobj, fdnames);
            
            % Assign the components : Electrolyte, NegativeElectrode, PositiveElectrode
            model.Electrolyte       = Electrolyte(paramobj.elyte);
            model.NegativeElectrode = Electrode(paramobj.ne);
            model.PositiveElectrode = Electrode(paramobj.pe);
            
        end
    
        function [problem, state] = getEquations(model, state0, state,dt, drivingForces, varargin)
            
            time = state0.time + dt;
            state = model.initStateAD(state);
            
            %% for now temperature and SOC are kept constant
            nc = model.G.cells.num;
            state.T   = model.T*ones(nc, 1);
            state.SOC = model.SOC*ones(nc, 1);
            
            % shortcuts
            battery = model;
            ne = 'NegativeElectrode';
            pe = 'PositiveElectrode';
            eac = 'ElectrodeActiveComponent';
            cc = 'CurrentCollector';
            elyte = 'Electrolyte';

            automaticAssembly;
            
            %% Set up the governing equations
            
            eqs={};
            
            %% We collect mass and charge conservation equations for the electrolyte and the electrodes

            names = {'Electrolyte', 'NegativeElectrode', 'PositiveElectrode'};
            
            for i = 1 : numel(names)
                eqs{end + 1} = model.getProp(state,{names{i}, 'massCons'});
                eqs{end + 1} = model.getProp(state,{names{i}, 'chargeCons'});
            end
            
            %% We collect charge conservation equations for the current collectors
            
            names = {'NegativeCurrentCollector', 'PositiveCurrentCollector'};
            for i = 1 : numel(names)
                eqs{end + 1} = model.getProp(state, {names{i}, 'chargeCons'});
            end
            
            %% We setup and add the control equation (fixed total current at PositiveCurrentCollector)
            
            src = drivingForces.src(time);
            coupterm = model.getCoupTerm('bc-PositiveCurrentCollector');
            faces = coupterm.couplingfaces;
            bcval = state.PositiveCurrentCollector.E;
            cond_pcc = model.PositiveCurrentCollector.EffectiveElectronicConductivity;
            [trans_pcc, cells] = model.PositiveCurrentCollector.operators.harmFaceBC(cond_pcc, faces);
            control = src - sum(trans_pcc.*(bcval - state.PositiveCurrentCollector.phi(cells)));
            
            eqs{end+1} = -control;


            
        end

        function state = updateT(model, state)
            names = {'NegativeElectrode', 'PositiveElectrode', 'Electrolyte'};
            for ind = 1 : numel(names);
                name = names{ind};
                nc = model.(name).G.cells.num;
                state.(name).T = state.T(1)*ones(nc, 1);
            end
        end
        
        function initstate = setupInitialState(model)
        % Setup initial state
        %
        % Abbreviations used in this function 
        % elyte : Electrolyte
        % ne    : NegativeElectrode
        % pe    : PositiveElectrode
        % ncc   : NegativeCurrentCollector
        % pcc   : PositiveCurrentCollector
            
            nc = model.G.cells.num;

            SOC = model.SOC;
            T   = model.T;
            
            initstate.T   =  T*ones(nc, 1);
            initstate.SOC =  SOC*ones(nc, 1);

            initstate = model.updateT(initstate);
            
            elyte = model.Electrolyte;
            ne    = model.NegativeElectrode;
            pe    = model.PositiveElectrode;
            ncc   = model.NegativeCurrentCollector;
            pcc   = model.PositiveCurrentCollector;
            ne_am = ne.ActiveMaterial;
            pe_am = pe.ActiveMaterial;

            %% setup initial NegativeElectrode state

            m = (1 ./ (ne_am.theta100 - ne_am.theta0));
            b = -m .* ne_am.theta0;
            theta = (SOC - b) ./ m;
            c = theta .* ne_am.Li.cmax;
            c = c*ones(ne.G.cells.num, 1);

            initstate.NegativeElectrode.ActiveMaterial.Li = c;
            initstate.NegativeElectrode.ActiveMaterial = ne_am.updateMaterialProperties(initstate.NegativeElectrode.ActiveMaterial);

            OCP = initstate.NegativeElectrode.ActiveMaterial.OCP;
            initstate.NegativeElectrode.ActiveMaterial.phi = OCP;

            %% setup initial PositiveElectrode state

            m = (1 ./ (pe_am.theta100 - pe_am.theta0));
            b = -m .* pe_am.theta0;
            theta = (SOC - b) ./ m;
            c = theta .* pe_am.Li.cmax;
            c = c*ones(pe.G.cells.num, 1);

            initstate.PositiveElectrode.ActiveMaterial.Li = c;
            initstate.PositiveElectrode.ActiveMaterial = pe_am.updateMaterialProperties(initstate.PositiveElectrode.ActiveMaterial);

            OCP = initstate.PositiveElectrode.ActiveMaterial.OCP;
            initstate.PositiveElectrode.ActiveMaterial.phi = OCP;

            %% setup initial Electrolyte state

            initstate.Electrolyte.phi = zeros(elyte.G.cells.num, 1);
            cs=cell(2,1);
            initstate.Electrolyte.cs = cs;
            initstate.Electrolyte.cs{1} = 1000*ones(elyte.G.cells.num, 1);

            %% setup initial Current collectors state
            
            OCP = initstate.NegativeElectrode.ActiveMaterial.OCP;
            OCP = OCP(1) .* ones(ncc.G.cells.num, 1);
            initstate.NegativeCurrentCollector.phi = OCP;

            OCP = initstate.PositiveElectrode.ActiveMaterial.OCP;
            OCP = OCP(1) .* ones(pcc.G.cells.num, 1);
            initstate.PositiveCurrentCollector.phi = OCP;
            
            initstate.PositiveCurrentCollector.E = OCP(1);
            
        end
        
        function state = setupElectrolyteCoupling(model, state)
        % Setup the electrolyte coupling by adding ion sources from the electrodes
        % shortcuts:
        % elyte : Electrolyte
        % ne : NegativeElectrode
        % pe : PositiveElectrode
            
            elyte = model.Electrolyte;
            ionSourceName = elyte.ionSourceName;
            coupterms = model.couplingTerms;
            
            phi = state.Electrolyte.phi;
            if isa(phi, 'ADI')
                adsample = getSampleAD(phi);
                adbackend = model.AutoDiffBackend;
                elyte_Li_source = adbackend.convertToAD(elyte_Li_source, adsample);
            end
            
            ne_R = state.NegativeElectrode.ElectrodeActiveComponent.ActiveMaterial.R;
            coupterm = getCoupTerm(couplingterms, 'NegativeElectrode-Electrolyte');
            elytecells = coupterm.couplingcells(:, 2);
            elyte_Li_source(elytecells) = ne_R;            
            
            pe_R = state.PositiveElectrode.ElectrodeActiveComponent.ActiveMaterial.R;
            coupterm = getCoupTerm(couplingterms, 'PositiveElectrode-Electrolyte');
            elytecells = coupterm.couplingcells(:, 2);
            elyte_Li_source(elytecells) = pe_R;
            
            state.Electrolyte.(ionSourceName) = elyte_Li_source;
        
        end
        
        
        function state = setupElectrodeCoupling(model, state)
        % Setup electrod coupling by updating the potential and concentration of the electrolyte in the active component of the
        % electrode. There, those quantities are considered as input and used to compute the reaction rate.
        %
        %
        % WARNING : at the moment, we do not pass the concentrations
        %
        % shortcuts:
        % elyte : Electrolyte
        % neac  : NegativeElectrode.ElectrodeActiveComponent 
        % peac  : PositiveElectrode.ElectrodeActiveComponent
            
            elyte = model.Electrolyte;
            neac = model.NegativeElectrode.ElectrodeActiveComponent;
            peac = model.PositiveElectrode.ElectrodeActiveComponent;
            
            phi_elyte = state.Electrolyte.phi;
            
            elyte_cells = zeros(model.G.cells.num, 1);
            elyte_cells(elyte.G.mappings.cellmap) = (1 : elyte.G.cells.num)';

            phi_elyte_neac = phi_elyte(elyte_cells(neac.G.mappings.cellmap));
            phi_elyte_peac = phi_elyte(elyte_cells(peac.G.mappings.cellmap));

            state.NegativeElectrode.ElectrodeActiveComponent.phiElectrolyte = phi_elyte_neac;
            state.PositiveElectrode.ElectrodeActiveComponent.phiElectrolyte = phi_elyte_peac;
            
        end
        
        
    end
    
end
