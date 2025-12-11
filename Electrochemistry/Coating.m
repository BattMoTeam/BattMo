classdef Coating < ElectronicComponent

    properties

        %% Sub-Models

        ActiveMaterial
        Binder
        ConductingAdditive

        % The two following models are instantiated only when activeMaterialModelSetup.composite is true, in this case,
        % ActiveMaterial model will remain empty. If activeMaterialModelSetup.composite is false, then the two models remains empty
        ActiveMaterial1
        ActiveMaterial2

        %% Input Parameters

        % Standard parameters
        effectiveDensity     % the mass density of the material (symbol: rho). Important : the density is computed with respect to total volume (including the empty pores)
        bruggemanCoefficient % the Bruggeman coefficient for effective transport in porous media (symbol: beta)
        
        activeMaterialModelSetup % Structure which describes the chose model, see schema in Utilities/JsonSchemas/Coating.schema.json. Here, we summarize
                              % - 'composite' : boolean (default is false)
                              % - 'SEImodel' : string with one of
                              %                 "none" (default)
                              %                 "Safari"
                              %                 "Bolay"
        
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
        compInds  % index of the sub models in the massFractions structure
        compnames % names of the components

        specificVolumes % One value per component (ActiveMaterial, Binder, ConductingAdditive) giving the specific volume (volume of 1kg of the component)
    end

    methods

        function model = Coating(inputparams)

            model = model@ElectronicComponent(inputparams);

            fdnames = {'effectiveDensity'               , ...
                       'bruggemanCoefficient'           , ...
                       'activeMaterialModelSetup'          , ...
                       'volumeFractions'                , ...
                       'volumeFraction'                 , ...
                       'thermalConductivity'            , ...
                       'specificHeatCapacity'           , ...
                       'effectiveThermalConductivity'   , ...
                       'effectiveVolumetricHeatCapacity', ...
                       'externalCouplingTerm'};

            model = dispatchParams(model, inputparams, fdnames);

            % Shortcuts
            sd  = 'SolidDiffusion';
            am  = 'ActiveMaterial';
            bd  = 'Binder';
            ad  = 'ConductingAdditive';
            sei = 'SolidElectrodeInterface';
            sr  = 'SideReaction';
            
            if model.activeMaterialModelSetup.composite
                am1 = 'ActiveMaterial1';
                am2 = 'ActiveMaterial2';
                compnames = {am1, am2, bd, ad};
            else 
                am = 'ActiveMaterial';
                compnames = {am, bd, ad};
            end

            model.compnames        = compnames;
            model.subModelNameList = compnames;

            %% We setup the volume fractions of each components

            % Compute the specific volumes of each component based on the density and mass fraction
            % If components or values are missing, the volume fraction is set to zero
            specificVolumes = zeros(numel(compnames), 1);
            for icomp = 1 : numel(compnames)
                compname = compnames{icomp};
                compInds.(compname) = icomp;
                if ~isempty(inputparams.(compname))
                    rho = inputparams.(compname).density;
                    mf  = inputparams.(compname).massFraction;
                    if ~isempty(rho) & ~isempty(mf)
                        specificVolumes(icomp) = mf/rho;
                    end
                end
            end

            model.compInds        = compInds;
            model.specificVolumes = specificVolumes;
            
            % We treat special cases for the specific volumes

            use_am_only = false;
            if model.activeMaterialModelSetup.composite
                if all(specificVolumes([compInds.(am1), compInds.(am2)]) == 0)
                    assert(~isempty(model.volumeFractions) && ~isempty(model.volumeFraction), ...
                           'Data in the subcomponents are missing. You can also provide volumeFractions and volumeFraction directly' )
                end
            else
                if all(specificVolumes == 0)
                    use_am_only = true;
                else
                    if specificVolumes(compInds.(am)) == 0
                        error('missing density and/or massFraction for the active material. The volume fraction cannot be computed ');
                    end
                end
            end

            % We normalize the volume fractions

            if isempty(model.volumeFractions)

                updateMassFractions = false;

                if use_am_only
                    % No data has been given, we assume that there is no binder and conducting additive
                    model.volumeFractions = zeros(numel(compnames), 1);
                    model.volumeFractions(compInds.(am)) = 1;
                    inputparams.(am).massFraction = 1;
                    inputparams.(bd).massFraction = 0;
                    inputparams.(ad).massFraction = 0;
                else                    
                    volumeFractions = zeros(numel(compnames), 1);
                    sumSpecificVolumes = sum(specificVolumes);
                    for icomp = 1 : numel(compnames)
                        volumeFractions(icomp) = specificVolumes(icomp)/sumSpecificVolumes;
                    end

                    model.volumeFractions = volumeFractions;
                end

            else

                updateMassFractions = true;

            end

            %% We setup the volume fraction of coating

            if isempty(model.volumeFraction)

                updateEffectiveDensity = false;

                % If the volumeFraction is given, we use it otherwise it is computed from the density and the specific
                % volumes of the components
                assert(~isempty(model.effectiveDensity), 'At this point we need an effective density in the model');
                model.volumeFraction = sumSpecificVolumes*model.effectiveDensity;

            else

                updateEffectiveDensity = true;

            end

            %% We setup the electronic conductivity

            if isempty(model.electronicConductivity)
                % if electronic conductivity is given, we use it otherwise we compute it as a volume average of the component
                kappa = 0;
                for icomp = 1 : numel(compnames)
                    compname = compnames{icomp};
                    if ~isempty(inputparams.(compname)) && ~isempty(inputparams.(compname).electronicConductivity)
                        kappa = kappa + model.volumeFractions(icomp)*inputparams.(compname).electronicConductivity;
                    end
                end
                assert(abs(kappa) > 0, 'The electronicConductivity must be provided for at least one component');
                model.electronicConductivity = kappa;
            end

            if isempty(model.effectiveElectronicConductivity)

                kappa = model.electronicConductivity;
                vf    = model.volumeFraction;
                bg    = model.bruggemanCoefficient;

                model.effectiveElectronicConductivity = kappa*vf^bg;

            end


            %% Setup the submodels

            np = inputparams.G.getNumberOfCells();

            if model.activeMaterialModelSetup.composite
                
                ams = {am1, am2};
                for iam = 1 : numel(ams)
                    amc = ams{iam};
                    inputparams.(amc).(sd).volumeFraction = model.volumeFraction*model.volumeFractions(model.compInds.(amc));
                    if strcmp(inputparams.(amc).diffusionModelType, 'full')
                        inputparams.(amc).(sd).np = np;
                    end
                    model.(amc) = ActiveMaterial(inputparams.(amc));
                end

            else
                
                switch model.activeMaterialModelSetup.SEImodel
                    
                  case {'none', 'Bolay'}

                    inputparams.(am).(sd).volumeFraction = model.volumeFraction*model.volumeFractions(model.compInds.(am));

                    switch inputparams.(am).diffusionModelType
                      case {'full', 'swelling'}
                        inputparams.(am).(sd).np = np;
                      case {'simple'}
                        % do nothing
                      otherwise
                        error('diffusion model type not recognized')
                    end
                    
                    model.ActiveMaterial = ActiveMaterial(inputparams.ActiveMaterial);
                    
                  case 'Safari'
                    
                    inputparams.(am).(sd).volumeFraction = model.volumeFraction*model.volumeFractions(model.compInds.(am));
                    inputparams.(am).(sd).np  = np;
                    inputparams.(am).(sei).np = np;
                    model.ActiveMaterial = SEIActiveMaterial(inputparams.ActiveMaterial);
                    
                  otherwise
                    
                    error('SEI model not recognized')
                    
                end
            end

            model.Binder             = Binder(inputparams.Binder);
            model.ConductingAdditive = ConductingAdditive(inputparams.ConductingAdditive);

            if updateMassFractions

                inputparams = model.updateMassFractions(inputparams);

            end

            if updateEffectiveDensity

                model = model.updateEffectiveDensity(inputparams);

            end


            %% We setup the thermal parameters

            if model.use_thermal

                %% We setup the thermal conductivities

                if isempty(model.thermalConductivity)
                    thermalConductivity = 0;
                    for icomp = 1 : numel(compnames)
                        compname = compnames{icomp};
                        if ~isempty(model.(compname).thermalConductivity)
                            thermalConductivity = thermalConductivity + model.volumeFractions(icomp)*inputparams.(compname).thermalConductivity;
                        end
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
                        if ~isempty(inputparams.(compname).specificHeatCapacity)
                            specificHeatCapacity = specificHeatCapacity + inputparams.(compname).massFraction*inputparams.(compname).specificHeatCapacity;
                        end
                    end
                    model.specificHeatCapacity = specificHeatCapacity;
                end

                if isempty(model.effectiveVolumetricHeatCapacity)
                    model.effectiveVolumetricHeatCapacity = model.volumeFraction*model.effectiveDensity*model.specificHeatCapacity;
                end

            end

        end

        function model = registerVarAndPropfuncNames(model)

            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@ElectronicComponent(model);

            itf = 'Interface';
            sd  = 'SolidDiffusion';
            sei = 'SolidElectrodeInterface';
            sr  = 'SideReaction';
            
            varnames = {'jCoupling', ...
                        'jExternal', ...
                        'SOC'};

            model = model.registerVarNames(varnames);

            % We declare SOC as an extra variable, as it is not used in assembly (otherwise it will be systematically
            % computed but not used)
            model = model.setAsExtraVarName('SOC');

            if ~model.activeMaterialModelSetup.composite

                am = 'ActiveMaterial';

                if strcmp(model.activeMaterialModelSetup.SEImodel, 'Bolay')
                    fn = @Coating.updateBolayEsource;
                    model = model.registerPropFunction({'eSource', fn, {{am, sd, 'Rvol'}, {am, itf, 'SEIflux'}}});
                else
                    fn = @Coating.updateEsource;
                    model = model.registerPropFunction({'eSource', fn, {{am, sd, 'Rvol'}}});
                end
                
                fn = @Coating.updatePhi;
                model = model.registerPropFunction({{am, itf, 'phiElectrode'}, fn, {'phi'}});

                fn = @Coating.dispatchTemperature;
                model = model.registerPropFunction({{am, 'T'}, fn, {'T'}});

                fn = @Coating.updateSOC;
                model = model.registerPropFunction({'SOC', fn, {{am, sd, 'cAverage'}}});

                if strcmp(model.activeMaterialModelSetup.SEImodel, 'Safari')
                    fn = @Coating.updateSideReactionPhi;
                    model = model.registerPropFunction({{am, sr, 'phiElectrode'}, fn, {'phi'}});
                end

            else
                  
                am1 = 'ActiveMaterial1';
                am2 = 'ActiveMaterial2';

                varnames = {{am1, 'SOC'}, ...
                            {am2, 'SOC'}};

                model = model.registerVarNames(varnames);
                model = model.setAsExtraVarNames(varnames);

                % We remove the dUdT variable (not used for non thermal simulation)
                ams = {am1, am2};
                for iam = 1 : numel(ams)
                    am = ams{iam};
                    if ~model.(am).(itf).includeEntropyChange
                        varname = {am, itf, 'dUdT'};
                        model = model.removeVarNames(varname);
                    end
                end

                fn = @Coating.updateCompositeEsource;
                inputnames = {{am1, sd, 'Rvol'}, ...
                              {am2, sd, 'Rvol'}};
                model = model.registerPropFunction({'eSource', fn, inputnames});

                fn = @Coating.updateCompositePhi;
                model = model.registerPropFunction({{am1, itf, 'phiElectrode'}, fn, {'phi'}});
                model = model.registerPropFunction({{am2, itf, 'phiElectrode'}, fn, {'phi'}});

                fn = @Coating.dispatchCompositeTemperature;
                model = model.registerPropFunction({{am1, 'T'}, fn, {'T'}});
                model = model.registerPropFunction({{am2, 'T'}, fn, {'T'}});

                fn = @Coating.updateCompositeSOC;
                inputnames = {{am1, sd, 'cAverage'}, ...
                              {am2, sd, 'cAverage'}};
                model = model.registerPropFunction({'SOC', fn, inputnames});
                model = model.registerPropFunction({{am1, 'SOC'}, fn, inputnames});
                model = model.registerPropFunction({{am2, 'SOC'}, fn, inputnames});

            end

            
            if model.use_thermal
                varnames = {'jFaceCoupling', ...
                            'jFaceExternal'};
                model = model.registerVarNames(varnames);

            end

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
                                                   
        end

        function model = setTPFVgeometry(model, tPFVgeometry)
        % tPFVgeometry should be instance of TwoPointFiniteVolumeGeometry

            model = setTPFVgeometry@ElectronicComponent(model, tPFVgeometry);

        end

        function inputparams = updateMassFractions(model, inputparams)

            compnames = model.compnames;

            for icomp = 1 : numel(compnames)
                compname = compnames{icomp};
                if ~isempty(inputparams.(compname).density)
                    massfractions(icomp) = inputparams.(compname).density*model.volumeFractions(icomp);
                else
                    massfractions(icomp) = 0;
                end
            end

            massfractions = massfractions/sum(massfractions);

            for icomp = 1 : numel(compnames)
                compname = compnames{icomp};
                inputparams.(compname).massFraction = massfractions(icomp);
            end

        end

        function jsonstruct = exportParams(model)

            jsonstruct = exportParams@ElectronicComponent(model);
            
            fdnames = {'effectiveDensity'            , ...     
                       'bruggemanCoefficient'        , ... 
                       'volumeFractions'             , ...                 
                       'volumeFraction'              , ...
                       'thermalConductivity'         , ...             
                       'specificHeatCapacity'        , ...            
                       'effectiveThermalConductivity', ...    
                       'effectiveVolumetricHeatCapacity' };
            
            for ifd = 1 : numel(fdnames)
                fdname = fdnames{ifd};
                jsonstruct.(fdname) = model.(fdname);
            end

        end
        
        
        function model = updateEffectiveDensity(model, inputparams)

            compnames = model.compnames;
            vf = model.volumeFraction;

            for icomp = 1 : numel(compnames)
                compname = compnames{icomp};
                if ~isempty(model.(compname).density)
                    massfractions(icomp) = inputparams.(compname).density*model.volumeFractions(icomp);
                else
                    massfractions(icomp) = 0;
                end
            end

            rho = sum(massfractions)*vf;

            model.effectiveDensity = rho;

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

        function state = updateEsource(model, state)

            am  = 'ActiveMaterial';
            sd  = 'SolidDiffusion';
            itf = 'Interface';

            F    = model.constants.F;
            n    = model.(am).(itf).numberOfElectronsTransferred;
            vols = model.G.getVolumes();

            state.eSource =  - n*F*vols.*state.(am).(sd).Rvol;

        end

        function state = updateBolayEsource(model, state)

            am  = 'ActiveMaterial';
            sd  = 'SolidDiffusion';
            itf = 'Interface';

            F    = model.constants.F;
            n    = model.(am).(itf).numberOfElectronsTransferred;
            vsa  = model.(am).(itf).volumetricSurfaceArea;
            
            vols = model.G.getVolumes();

            Rvol    = state.(am).(sd).Rvol;
            seiflux = state.(am).(itf).SEIflux;
            
            state.eSource =  F*vols.*( -n*Rvol + vsa*seiflux );

        end

        
       function state = updateCompositeEsource(model, state)

            am1 = 'ActiveMaterial1';
            am2 = 'ActiveMaterial2';
            itf = 'Interface';
            sd  = 'SolidDiffusion';

            F    = model.constants.F;
            vols = model.G.getVolumes();

            ams = {am1, am2};

            Rvol = 0;

            for iam = 1 : numel(ams)

                amc = ams{iam};

                n    = model.(amc).(itf).numberOfElectronsTransferred;
                Rvol = Rvol + n*F*state.(amc).(sd).Rvol;

            end

            state.eSource = - vols.*Rvol;

        end


        function state = updateCompositePhi(model, state)

            am1 = 'ActiveMaterial1';
            am2 = 'ActiveMaterial2';
            itf = 'Interface';

            ams = {am1, am2};

            for iam = 1 : numel(ams)

                amc = ams{iam};
                state.(amc).(itf).phiElectrode = state.phi;

            end

        end

        function state = updatePhi(model, state)

            am  = 'ActiveMaterial';
            itf = 'Interface';

            state.(am).(itf).phiElectrode = state.phi;

        end

        function state = updateSideReactionPhi(model, state)

            am  = 'ActiveMaterial';
            sr = 'SideReaction';

            state.(am).(sr).phiElectrode = state.phi;
            
        end
        
        function state = dispatchCompositeTemperature(model, state)

            am1 = 'ActiveMaterial1';
            am2 = 'ActiveMaterial2';

            ams = {am1, am2};

            for iam = 1 : numel(ams)

                amc = ams{iam};
                state.(amc).T = state.T;

            end

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

            vf       = model.volumeFraction;
            am_frac  = model.volumeFractions(model.compInds.(am));
            vols     = model.G.getVolumes();
            cmax     = model.(am).(itf).saturationConcentration;
            theta100 = model.(am).(itf).guestStoichiometry100;
            theta0   = model.(am).(itf).guestStoichiometry0;

            c = state.(am).(sd).cAverage;
            
            %% We do not use the gueststochiometry value to compute the State of Charge
            
            theta = c/cmax;
            % m     = (1 ./ (theta100 - theta0));
            % b     = -m .* theta0;
            % SOC   = theta*m + b;
            SOC = theta;
            vol = am_frac*vf.*vols;

            SOC = sum(SOC.*vol)/sum(vol);

            state.SOC = SOC;

        end

        function state = updateCompositeSOC(model, state)

            am1 = 'ActiveMaterial1';
            am2 = 'ActiveMaterial2';
            itf = 'Interface';
            sd  = 'SolidDiffusion';

            vf    = model.volumeFraction;
            vols = model.G.getVolumes();

            ams = {am1, am2};

            for iam = 1 : numel(ams)

                amc = ams{iam};

                am_frac  = model.volumeFractions(model.compInds.(amc));
                cmax     = model.(amc).(itf).saturationConcentration;
                theta100 = model.(amc).(itf).guestStoichiometry100;
                theta0   = model.(amc).(itf).guestStoichiometry0;

                c = state.(amc).(sd).cAverage;

                vol = am_frac*vf.*vols;

                %% We do not use the gueststochiometry value to compute the State of Charge
                
                molvals(iam)    = sum(c.*vol);
                % molval0s(iam)   = theta0*cmax*sum(vol);
                % molval100s(iam) = theta100*cmax*sum(vol);

                molval0s(iam)   = 0;
                molval100s(iam) = cmax*sum(vol);
                
                state.(amc).SOC = (molvals(iam) - molval0s(iam))/(molval100s(iam) - molval0s(iam));

            end

            molval    = sum(molvals);
            molval0   = sum(molval0s);
            molval100 = sum(molval100s);

            state.SOC = (molval - molval0)/(molval100 - molval0);

        end


    end
end



%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
