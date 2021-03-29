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
        
        function model = Battery2(params)
            
            model = PhysicalModel([]);
            
            %% Setup the model using the input parameters
            model.T     = params.T;
            model.SOC   = params.SOC;
            model.J     = params.J;
            model.Ucut  = params.Ucut;
            
            % Assign the grid
            model.G = params.G;
            
            % Assign the autodiff backend
            model.G = params.AutoDiffBackend;            
            
            % Assign the components : Electrolyte, NegativeElectrode, PositiveElectrode
            model.Electrolyte       = params.Electrolyte;
            model.NegativeElectrode = params.NegativeElectrode;
            model.PositiveElectrode = params.PositiveElectrode;
            
            % Assign the coupling term structure
            model.couplingTerms = params.couplingTerms;
            
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
