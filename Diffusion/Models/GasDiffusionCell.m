classdef GasDiffusionCell < MaxwellStefanGasDiffusion

    properties

        externalCouplingTerms % Cell array with external coupling terms

        Control
        
    end
    
    methods

        function model = GasDiffusionCell(inputparams)

            model = model@MaxwellStefanGasDiffusion(inputparams);
    
            fdnames = {'externalCouplingTerms'};
            model = dispatchParams(model, inputparams, fdnames);

            model.Control = GasDiffusionCellControl(inputparams.Control);

            model.subModelNameList{end + 1} = 'Control';

            model = model.setupControlMappings();
            
        end

        function model = setupControlMappings(model)

            nctrl = model.Control.nControls;

            faces = [];
            ctrl  = [];
            
            for ictrl = 1 : nctrl

                coupbcfaces = model.externalCouplingTerms{ictrl}.couplingfaces;

                faces = vertcat(faces, coupbcfaces);
                ctrl  = vertcat(ctrl, repmat(ictrl, numel(coupbcfaces), 1));
                
            end

            bcfacectrltbl.faces = faces;
            bcfacectrltbl.ctrl  = ctrl;
            bcfacectrltbl = IndexArray(bcfacectrltbl);

            
            ctrltbl.ctrl = (1 : nctrl)';
            ctrltbl = IndexArray(ctrltbl);

            map = TensorMap();
            map.fromTbl  = ctrltbl;
            map.toTbl    = bcfacectrltbl;
            map.mergefds = {'ctrl'};
            map = map.setup();

            ctrlToBc = map.getMatrix();
            bcToCtrl = ctrlToBc';

            model.boundaryFaces = bcfacectrltbl.get('faces');
            model.Control.mappings.ctrlToBc = ctrlToBc;
            model.Control.mappings.bcToCtrl = bcToCtrl;
            
            model = model.setupMappings();
            
        end
        
        function model = registerVarAndPropfuncNames(model)

            ncomp = model.numberOfComponents;
            
            model = registerVarAndPropfuncNames@MaxwellStefanGasDiffusion(model);

            model = model.setAsStaticVarName(VarName({}, 'sources', ncomp));
            
            fn = @GasDiffusionCell.updateControlValue;
            model = model.registerPropFunction({{'Control', 'value'}, fn, {{'Control', 'type'}}});

            fn = @GasDiffusionCell.updateControlFluxBoundaryEquation;
            inputvarnames = {VarName({'Boundary'}, 'diffusionFluxes', ncomp), ...
                             VarName({'Control'}, 'flux')};
            outputvarnames = VarName({'Control'}, 'fluxBoundaryEquations', ncomp);
            model = model.registerPropFunction({outputvarnames, fn, inputvarnames});

            for icomp = 1 : ncomp
                fn = @GasDiffusionCell.updateBoundaryMassFractions;
                inputvarnames = {VarName({'Control'}, 'massFractions', ncomp, icomp)};
                outputvarnames = VarName({'Boundary'}, 'massFractions', ncomp, icomp);
                model = model.registerPropFunction({outputvarnames, fn, inputvarnames});
            end
            
            fn = @GasDiffusionCell.updateControlBoundaryPressureEquation;
            inputvarnames = {{'Control', 'pressure'}, ...
                             {'Boundary', 'pressure'}};
            outputvarname = {'Control', 'pressureBoundaryEquation'};
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});

        end

        function initstate = setupInitialState(model, jsonstruct)

            nc    = model.G.getNumberOfCells();
            nbf   = numel(model.boundaryFaces);
            nctrl = model.Control.nControls;
            
            ncomp = model.numberOfComponents;
            mws   = model.molarWeights;
            
            ps = jsonstruct.initialPressures;

            % Setup initial pressure
            initstate.pressure = sum(ps)*ones(nc, 1);
            initstate.Boundary.pressure = sum(ps)*ones(nbf, 1);
            
            % compute mass fractions
            mfs = 1/sum(mws.*ps)*(mws.*ps);

            % Setup mass fraction
            for icomp = 1 : ncomp
                initstate.massFractions{icomp}            = mfs(icomp)*ones(nc, 1);
                initstate.Boundary.massFractions{icomp}   = mfs(icomp)*ones(nbf, 1);
                initstate.Boundary.diffusionForces{icomp} = zeros(nbf, 1);
            end

            initstate.Control.pressure = sum(ps)*ones(nctrl, 1);
            initstate.Control.flux     = zeros(nctrl, 1);

            ctrlelts = model.Control.controlElements;

            type  = cellfun(@(ctrlelt) ctrlelt.type, ctrlelts);
            value = cellfun(@(ctrlelt) ctrlelt.type, ctrlelts);

            initstate.Control.pressure(type == 1) = value(type == 1);
            initstate.Control.flux(type == 2)     = value(type == 2);
            
        end
        
        function state = updateControlValue(model, state)

            ctrlelts = model.Control.controlElements

            state.Control.values = cellfun(@(elt) elt.value, ctrlelts);
            
        end

        function state = updateControlBoundaryPressureEquation(model, state)

            ctrlToBc = model.Control.mappings;

            pbc = state.Boundary.pressure;
            pc  = state.Control.pressure;
            
            state.Control.pressureBoundaryEquation = ctrlToBc*pc - pbc;
            
        end

        function state = updateBoundaryMassFractions(model, state)

            ncomp = model.numberOfComponents;

            ctrlToBc = model.Control.mappings.ctrlToBc;

            for icomp = 1 : ncomp

                state.Boundary.massFractions{icomp} = ctrlToBc*state.Control.massFractions{icomp};

            end
            
        end
        
        
    end

end
