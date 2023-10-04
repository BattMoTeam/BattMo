classdef Coating < ElectronicComponent

    properties

        %% Sub-Models

        ActiveMaterial
        Binder
        ConductingAdditive

        %% Input Parameters

        % Standard parameters
        density              % the mass density of the material (symbol: rho). Important : the density is computed with respect to total volume (including the empty pores)
        bruggemanCoefficient % the Bruggeman coefficient for effective transport in porous media (symbol: beta)
        activematerial_type  % 'default' (only one particle type) or 'composite' (two different particles)
        
        % Advanced parameters (used if given, otherwise computed)
        volumeFractions                 % mass fractions of each components (if not given computed subcomponent and density)
        volumeFraction
        thermalConductivity             % (if not given computed from the subcomponents)
        specificHeatCapacity            % (if not given computed from the subcomponents)
        effectiveThermalConductivity    % (account for volume fraction, if not given, computed from thermalConductivity)
        effectiveVolumetricHeatCapacity % (account for volume fraction and density, if not given, computed from specificHeatCapacity)

        % external coupling parameters
        externalCouplingTerm % structure to describe external coupling (used in absence of current collector)
        
        %% Computed parameters at model setup
        compInds            % index of the sub models in the massFractions structure
        
    end

    methods
        
        function model = Coating(paramobj)

            model = model@ElectronicComponent(paramobj);
            
            fdnames = {'density'                        , ...
                       'bruggemanCoefficient'           , ...
                       'activematerial_type'            , ...
                       'volumeFractions'                , ...
                       'volumeFraction'                 , ...
                       'thermalConductivity'            , ...
                       'specificHeatCapacity'           , ...
                       'effectiveThermalConductivity'   , ...
                       'effectiveVolumetricHeatCapacity', ...
                       'externalCouplingTerm'};

            model = dispatchParams(model, paramobj, fdnames);

            am = 'ActiveMaterial';
            sd = 'SolidDiffusion';
            am = 'ActiveMaterial';
            bd = 'Binder';
            ad = 'ConductingAdditive';

            compnames = {am, bd, ad};

            %% We setup the volume fractions of each components
            
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

            if isempty(model.volumeFractions)

                sumSpecificVolumes = sum(specificVolumes);
                for icomp = 1 : numel(compnames)
                    volumeFractions(icomp) = specificVolumes(icomp)/sumSpecificVolumes;
                end

                model.volumeFractions = volumeFractions;
                
            end

            %% We setup the volume fraction of coating
            
            if isempty(model.volumeFraction)
                
                % If the volumeFraction is given, we use it otherwise it is computed from the density and the specific
                % volumes of the components
                model.volumeFraction = sumSpecificVolumes*model.density;
                
            end

            %% We setup the electronic conductivity

            if isempty(model.electronicConductivity)
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

            if isempty(model.effectiveElectronicConductivity)

                kappa = model.electronicConductivity;
                vf    = model.volumeFraction;
                bg    = model.bruggemanCoefficient;
                
                model.effectiveElectronicConductivity = kappa*vf^bg;
                
            end

            %% We setup the thermal parameters
            
            if model.use_thermal
                
                %% We setup the thermal conductivities
                
                if isempty(model.thermalConductivity)
                    bg = model.bruggemanCoefficient;
                    thermalConductivity = 0;
                    for icomp = 1 : numel(compnames)
                        compname = compnames{icomp};
                        thermalConductivity = thermalConductivity + (model.volumeFractions(icomp))^bg*paramobj.(compname).thermalConductivity;
                    end
                    model.thermalConductivity = thermalConductivity;
                end

                if isempty(model.effectiveThermalConductivity)
                    bg = model.bruggemanCoefficient;
                    model.effectiveThermalConductivity = (model.volumeFraction).^bg.*model.thermalConductivity;
                end

                %% We setup the thermal capacities

                if isempty(model.specificHeatCapacity)
                    specificHeatCapacity = 0;
                    for icomp = 1 : numel(compnames)
                        compname = compnames{icomp};
                        specificHeatCapacity = specificHeatCapacity + paramobj.(compname).massFraction*paramobj.(compname).specificHeatCapacity;
                    end
                    model.specificHeatCapacity = specificHeatCapacity;
                end

                if isempty(model.effectiveVolumetricHeatCapacity)
                    model.effectiveVolumetricHeatCapacity = model.density*specificHeatCapacity;
                end

            end
            
            %% Setup the submodels
            
            np = paramobj.G.cells.num;
            switch paramobj.activematerial_type
              case 'default'
                if strcmp(paramobj.(am).diffusionModelType, 'full')
                    paramobj.(am).(sd).volumeFraction = model.volumeFraction*model.volumeFractions(model.compInds.(am));
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

        function model = registerVarAndPropfuncNames(model)

            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@ElectronicComponent(model);

            am  = 'ActiveMaterial';
            itf = 'Interface';
            sd  = 'SolidDiffusion';
            
            varnames = {'jCoupling', ...
                        'jExternal', ...
                        'SOC'};

            model = model.registerVarNames(varnames);            

            if model.use_thermal
                varnames = {'jFaceCoupling', ...
                            'jFaceExternal'};
                model = model.registerVarNames(varnames);

            end
            
            fn = @Coating.updateCurrentSource;
            model = model.registerPropFunction({'eSource', fn, {{am, sd, 'Rvol'}}});
            
            fn = @Coating.updatePhi;
            model = model.registerPropFunction({{am, itf, 'phiElectrode'}, fn, {'phi'}});
            
            fn = @Coating.dispatchTemperature;
            model = model.registerPropFunction({{am, 'T'}, fn, {'T'}});

            fn = @Coating.updatejBcSource;
            model = model.registerPropFunction({'jBcSource', fn, {'jCoupling', 'jExternal'}});

            if model.use_thermal
                fn = @Coating.updatejFaceBc;
                model = model.registerPropFunction({'jFaceBc', fn, {'jFaceCoupling', 'jFaceExternal'}});
            end
            
            fn = @Coating.updatejExternal;
            model = model.registerPropFunction({'jExternal', fn, {}});
            if model.use_thermal
                model = model.registerPropFunction({'jFaceExternal', fn, {}});
            end
            
            fn = @Coating.updatejCoupling;
            model = model.registerPropFunction({'jCoupling', fn, {}});
            if model.use_thermal
                model = model.registerPropFunction({'jFaceCoupling', fn, {}});
            end

            fn = @Coating.updateSOC;
            model = model.registerPropFunction({'SOC', fn, {{am, sd, 'cAverage'}}});
            
            % We declare SOC as an extra variable, as it is not used in assembly (otherwise it will be systematically
            % computed but not used)
            model = model.setAsExtraVarName('SOC');

        end
        
        
        function state = updatejBcSource(model, state)
            
            state.jBcSource = state.jCoupling + state.jExternal;
            
        end
        
        function state = updatejFaceBc(model, state)
            
            state.jFaceBc = state.jFaceCoupling + state.jFaceExternal;
            
        end
        
        function state = updatejExternal(model, state)
            
            state.jExternal     = 0;
            state.jFaceExternal = 0;
            
        end

        function state = updatejCoupling(model, state)
            
            state.jCoupling     = 0;
            state.jFaceCoupling = 0;
            
        end
        
        function state = updateCurrentSource(model, state)
            
            am  = 'ActiveMaterial';
            itf = 'Interface';
            sd  = 'SolidDiffusion';
            
            F    = model.(am).(itf).constants.F;
            vols = model.G.cells.volumes;
            n    = model.(am).(itf).numberOfElectronsTransferred;

            Rvol = state.(am).(sd).Rvol;
            
            state.eSource = - vols.*Rvol*n*F; % C/s
            
        end
        
        function state = updatePhi(model, state)
            
            am  = 'ActiveMaterial';
            itf = 'Interface';
            
            state.(am).(itf).phiElectrode = state.phi;
            
        end

        function state = dispatchTemperature(model, state)

            am  = 'ActiveMaterial';

            state.(am).T = state.T;
            
        end

        
        function state = updateSOC(model, state)

            % shortcut
            am  = 'ActiveMaterial';
            itf = 'Interface';
            sd  = 'SolidDiffusion';

            vf       = model.(am).volumeFraction;
            am_frac  = model.(am).volumeFractions(model.compInds.(am));
            vols     = model.G.cells.volumes;
            cmax     = model.(itf).(itf).saturationConcentration;
            theta100 = model.(itf).(itf).guestStoichiometry100;
            theta0   = model.(itf).(itf).guestStoichiometry0;
            
            c = state.(am).(sd).cAverage;

            theta = c/cmax;
            m     = (1 ./ (theta100 - theta0));
            b     = -m .* theta0;
            SOC   = theta*m + b;
            vol   = am_frac*vf.*vols;
            
            SOC = sum(SOC.*vol)/sum(vol);

            state.SOC = SOC;
            
        end

    end
end

