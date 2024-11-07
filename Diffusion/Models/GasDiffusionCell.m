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

            mappings = struct('ctrlToBc', ctrlToBc, ...
                              'bcToCtrl', bcToCtrl);
            
            model.boundaryFaces    = bcfacectrltbl.get('faces');
            model.Control.mappings = mappings;
            
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

            fn = @GasDiffusionCell.updateBoundaryMassFractions;
            inputvarnames = {VarName({'Control'}, 'massFractions', ncomp)};
            outputvarnames = VarName({'Boundary'}, 'massFractions', ncomp);
            model = model.registerPropFunction({outputvarnames, fn, inputvarnames});
            
            fn = @GasDiffusionCell.updateControlBoundaryPressureEquation;
            inputvarnames = {{'Control', 'pressure'}, ...
                             {'Boundary', 'pressure'}};
            outputvarname = {'Control', 'pressureBoundaryEquation'};
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});

        end
        
    end

end
