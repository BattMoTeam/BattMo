classdef PorousTransportLayerBoundary < BaseModel
    
    properties
        
        compInd   % mapping structure for component indices
        phaseInd  % mapping structure for phase indices
        liquidInd % mapping structure for component indices
        gasInd    % mapping structure for component indices
        
    end
    
    methods
        
        function model = PorousTransportLayerBoundary()

            model = model@BaseModel();
            
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

            % Masses for each component of gas (per total volume, as used in mass of conservation law) in [kg m^-3]
            densities = VarName({}, 'gasDensities', ngas);
            varnames{end + 1} = densities;

            % Liquid density (Mass of liquid per volume of liquid) in [kg m^-3]
            varnames{end + 1} = 'liqrho';

            % concentration of OH in [mol m^-3]
            varnames{end + 1} = 'cOH';

            model = model.registerVarNames(varnames);
            
            model = model.removeVarName(VarName({}, 'phasePressures', nph, phaseInd.solid));
        end
        
    end


end

