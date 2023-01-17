classdef PorousTransportLayerBoundary < BaseModel
    
    properties

        compInd   % mapping structure for component indices
        phaseInd  % mapping structure for phase indices
        mobPhaseInd % mapping structure for mobile phase indices
        liquidInd % mapping structure for component indices
        gasInd    % mapping structure for component indices

        controlValues % Structure with following fields
                      % - cOH
                      % - liqrho
                      % - phasePressures
                      % - gasDensities
                      % - phi
        
    end
    
    methods
        
        function model = PorousTransportLayerBoundary(paramobj)
            
            model = model@BaseModel();

            fdnames = {'compInd'    , ...
                       'phaseInd'   , ...
                       'liquidInd'  , ...
                       'mobPhaseInd', ...
                       'gasInd'};
            model = dispatchParams(model, paramobj, fdnames);

            
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

            % Boundary equations
            bcEquations = VarName({}, 'bcEquations', 2 + ngas);
            varnames{end + 1} = bcEquations;
            
            % Boundary control equations
            bcControlEquations = VarName({}, 'bcControlEquations', nmobph + 1);
            varnames{end + 1} = bcControlEquations;

            model = model.registerVarNames(varnames);

            model = model.removeVarName(VarName({}, 'phasePressures', nph, phaseInd.solid));

            cOH = VarName({}, 'concentrations', nliquid, liquidInd.OH);
            fn = @() PorousTransportLayerBoundary.updateConcentrations;
            inputnames = {'liqrho', ...
                          cOH};
            ind = setdiff([1 : nliquid]', liquidInd.OH);
            model = model.registerPropFunction({VarName({}, 'concentrations', nliquid, ind), fn, inputnames});
            
            fn = @() PorousTransportLayerBoundary.setupBcControlEquations;
            inputnames = {VarName({}, 'phasePressures', nph, model.mobPhaseInd.phaseMap), ...
                          'phi'};
            model = model.registerPropFunction({bcControlEquations, fn, inputnames});
            
        end

        function state = updateConcentrations(model, state)

            state = PorousTransportLayer.updateConcentrations(model, state);
            
        end

        function state = setupBcControlEquations(model, state)

            nmobph = model.mobPhaseInd.nmobphase;
            
            for imobph = 1 : nmobph
                iph = model.mobPhaseInd.phaseMap(imobph);
                eqs{imobph} = state.phasePressures{iph} - model.controlValues.phasePressure{imobph};
            end

            eqs{nmobph + 1} = state.phi - model.controlValues.phi;

            state.bcControlEquations = eqs;
            
        end
        
        
    end


end

