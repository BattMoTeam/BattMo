classdef Coating < ElectronicComponent

    properties

        %% Sub-Models

        ActiveMaterial
        Binder
        ConductingAdditive

        %% Input Parameters

        % Standard parameters
        density              % the mass density of the material (symbol: rho)
        bruggemanCoefficient % the Bruggeman coefficient for effective transport in porous media (symbol: beta)
        activematerial_type % 'default' (only one particle type) or 'composite' (two different particles)
        
        % Advanced parameters (used if given, otherwise computed)
        effectiveThermalConductivity
        volumeFraction
        
        %% Computed parameters at model setup
        volumeFractions     % mass fractions of each components
        compInds            % index of the sub models in the massFractions structure
        
    end

    methods
        
        function model = Coating(paramobj)

            model = model@ElectronicComponent(paramobj);
            
            fdnames = {'density'                     , ...
                       'bruggemanCoefficient'        , ...
                       'activematerial_type'         , ...
                       'effectiveThermalConductivity', ...
                       'volumeFractions'             , ...
                       'volumeFraction'};

            model = dispatchParams(model, paramobj, fdnames);

            am = 'ActiveMaterial';
            sd = 'SolidDiffusion';
            am = 'ActiveMaterial';
            bd = 'Binder';
            ad = 'ConductingAdditive';

            compnames = {am, bd, ad};

            %% We compute the volume fractions of each components
            
            % Compute the specific volumes of each component based on the density and mass fraction
            % If components or values are missing, the volume fraction is set to zero
            for icomp = 1 : numel(compnames)
                compname = compnames{icomp};
                compInds.(compname) = icomp;
                if ~isempty(paramobj.(compname))
                    rho = paramobj.(compname).density;
                    mf  = paramobj.(compname).massFraction;
                    if ~isempty(rho) & ~isempty(mf)
                        specificVolumes(icomp) = mf/rho;
                    else
                        specificVolumes(icomp) = 0;
                    end
                else
                    specificVolumes(icomp) = 0;                    
                end
            end

            model.compInds = compInds;
            
            % We treat special cases
            
            if all(specificVolumes == 0)
                % No data has been given, we assume that there is no binder and conducting additive
                specificVolumes(compInds.(am)) = 1;
            else
                if specificVolumes(compInds.(am)) == 0
                    error('missing density and/or massFraction for the active material. The volume fraction cannot be computed ');
                end
            end

            % We normalize the volume fractions
            
            sumSpecificVolumes = sum(specificVolumes);
            for icomp = 1 : numel(compnames)
                volumeFractions(icomp) = specificVolumes(icomp)/sumSpecificVolumes;
            end

            model.volumeFractions = volumeFractions;

            %% We compute volume fraction of coating
            
            if isempty(paramobj.volumeFraction)
                % If the volumeFraction is given, we use it otherwise it is computed from the density and the specific
                % volumes of the components
                model.volumeFraction = sumSpecificVolumes*paramobj.density;
            end

            %% We compute the electronic conductivity

            if isempty(paramobj.electronicConductivity)
                % if electronic conductivity is given, we use it otherwise we compute it as a volume average of the component
                kappa = 0;
                for icomp = 1 : numel(compnames)
                    compname = compnames{icomp};
                    if ~isempty(paramobj.(compname))
                        assert(~isempty(paramobj.(compname)), 'missing electronic conductivity');
                        kappa = kappa + model.volumeFractions(icomp)*paramobj.(compname).electronicConductivity;
                    end
                end

                model.electronicConductivity = kappa;
            end

            if isempty(paramobj.effectiveElectronicConductivity)
                kappa = paramobj.electronicConductivity;
                vf    = paramobj.volumeFraction;
                bcoef = paramobj.bruggemanCoefficient; 
                model.effectiveElectronicConductivity = kappa*vf^bcoef;
            end

            %% Setup the submodels
            
            np = paramobj.G.cells.num;
            switch paramobj.activematerial_type
              case 'default'
                if strcmp(paramobj.(am).diffusionModelType, 'full')
                    paramobj.(am).(sd).np = np;
                end
                model.ActiveMaterial = ActiveMaterial(paramobj.ActiveMaterial);
              case 'composite'
                model.ActiveMaterial = CompositeActiveMaterial(paramobj.ActiveMaterial);
              otherwise
                error('activematerial_type not recognized');
            end
            
            model.Binder             = Binder(paramobj.Binder);
            model.ConductingAdditive = ConductingAdditive(paramobj.ConductingAdditive);
            
        end
        
    end
end

