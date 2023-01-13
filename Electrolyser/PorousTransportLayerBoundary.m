classdef PorousTransportLayerBoundary < BaseModel
    
    properties

        constants

        MWs % molecular weight for each of the gas components
        V0s % partial molar volumes
        
        compInd   % mapping structure for component indices
        phaseInd  % mapping structure for phase indices
        liquidInd % mapping structure for component indices
        gasInd    % mapping structure for component indices
        
    end
    
    methods
        
        function model = PorousTransportLayerBoundary(paramobj)
            
            model = model@BaseModel();

            model.MWs = paramobj.MWs;
            
            compInd.H2Oliquid = 1;
            compInd.H2Ogas    = 2;
            compInd.OH        = 3;
            compInd.K         = 4;
            compInd.activeGas = 5;
            % compInd.(H2 or O2) = 5 % should be instantiated by derived class see HydrogenPorousTransportLayer.m and OxygenPorousTransportLayer.m
            compInd.ncomp     = 5;
            compInd.liquid    = [compInd.H2Oliquid; compInd.OH; compInd.K];
            compInd.gas       = [compInd.H2Ogas; compInd.activeGas];
            
            phaseInd.liquid = 1;
            phaseInd.gas    = 2;
            phaseInd.solid  = 3;
            phaseInd.mobile = [1; 2];
            phaseInd.nphase = 3;
            
            % compInd.phaseMap(compInd.H2Oliquid)  = phaseInd.liquid;
            compInd.phaseMap  = [1; 2; 1; 1; 2]; % first component (H2Oliquid) is in phase indexed by 1 (liquid phase), and so on
            
            liquidInd.H2Oliquid = 1;
            liquidInd.OH = 2;
            liquidInd.K  = 3;
            liquidInd.ncomp  = 3;
            liqudInd.compMap = [1; 3; 4];
            
            gasInd.H2Ogas    = 1;
            gasInd.activeGas = 2;
            gasInd.ncomp     = 2;
            gasInd.compMap   = [2; 5];

            model.compInd = compInd;
            model.phaseInd = phaseInd;
            model.liquidInd = liquidInd;            
            model.gasInd = gasInd;
            
            model.constants = PhysicalConstants();
            
        end

        function model = registerVarAndPropfuncNames(model)
            
        %% Declaration of the Dynamical Variables and Function of the model
        % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@BaseModel(model);

            phaseInd  = model.phaseInd;
            liquidInd = model.liquidInd;
            gasInd    = model.gasInd;
            compInd   = model.compInd;

            ncomp   = compInd.ncomp;
            ngas    = gasInd.ncomp;
            nph     = phaseInd.nphase;
            nliquid = liquidInd.ncomp;
            nmobph  = numel(phaseInd.mobile);

            varnames = {};

            % Phase pressures in [Pa]
            phasePressures = VarName({}, 'phasePressures', nph);
            varnames{end + 1} = phasePressures;

            % Gas densities in mass per gas volume (one value for each gas component)
            gasDensities = VarName({}, 'gasDensities', ngas);
            varnames{end + 1} = gasDensities;
            
            % Liquid density (Mass of liquid per volume of liquid) in [kg m^-3]
            varnames{end + 1} = 'liqrho';

            % Concentrations in the liquid in [mol m^-3]
            concentrations = VarName({}, 'concentrations', nliquid);
            varnames{end + 1} = concentrations;

            % electric potential
            varnames{end + 1} = 'phi';

            % Temperature
            varnames{end + 1} = 'T';
            
            % Based on ideal gas
            varnames{end + 1} = 'gasStateEquation';

            % State equation for the liquid (for the moment, we assume incompressibility)
            varnames{end + 1} = 'liquidStateEquation';

            model = model.registerVarNames(varnames);

            fn = @PorousTransportLayerBoundary.updateGasStateEquation;
            inputnames = {'T'         , ...
                          gasDensities, ...
                          VarName({}, 'phasePressures', nph, phaseInd.gas)};
            model = model.registerPropFunction({'gasStateEquation', fn, inputnames});

            model = model.removeVarName(VarName({}, 'phasePressures', nph, phaseInd.solid));
            
            fn = @() PorousTransportLayerBoundary.liquidStateEquation;
            inputnames = {concentrations};
            model = model.registerPropFunction({'liquidStateEquation', fn, inputnames});

            fn = @() PorousTransportLayerBoundary.updateConcentrations;
            inputnames = {'liqrho', ...
                          VarName({}, 'concentrations', nliquid, liquidInd.OH)};
            ind = setdiff([1 : nliquid]', liquidInd.OH);
            model = model.registerPropFunction({VarName({}, 'concentrations', nliquid, ind), fn, inputnames});

        end

        function updateGasStateEquation(model, state)
        % We use ideal gas law

            R = model.constants.R;

            T = state.T;
            eq = state.phasePressures{model.phaseInd.gas};
            for igas = 1 : model.gasInd.ngas
                eq = eq - R*T.*state.gasDensities{igas}./MWs{igas};
            end
            
        end

        function state = liquidStateEquation(model, state)

            state = PorousTransportLayer.liquidStateEquation(model, state);
            
        end

        function state = updateConcentrations(model, state)

            state = PorousTransportLayer.updateConcentrations(model, state);
            
        end
        
        
    end


end

