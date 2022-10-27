classdef ElectronicComponent < BaseModel
%
% The ElectronicComponent class is model to assemble charge conservation equation
% 
    properties
        
        EffectiveElectricalConductivity % Effective electrical conductivity
        constants

        use_thermal
        
    end

    methods
        
        function model = ElectronicComponent(paramobj)
        % Here, :code:`paramobj` is instance of :class:`ElectronicComponentInputParams <Electrochemistry.ElectronicComponentInputParams>`

            model = model@BaseModel();
            
            % OBS : All the models should have same backend (this is not assigned automaticallly for the moment)
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks', false);
            
            fdnames = {'G'                              , ...
                       'EffectiveElectricalConductivity', ...
                       'use_thermal'};
            
            model = dispatchParams(model, paramobj, fdnames);
            
            % setup discrete differential operators
            docellflux = false;
            if model.use_thermal
                docellflux = true;
            end
            
            model.operators = localSetupOperators(model.G, 'assembleCellFluxOperator', docellflux);
            
            model.constants = PhysicalConstants();

        end

        function model = registerVarAndPropfuncNames(model)
            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@BaseModel(model);
            
            varnames = {'T'           , ...
                        'phi'         , ...
                        'jBcSource'   , ...
                        'eSource'     , ...
                        'conductivity', ...
                        'j'           , ...
                        'chargeCons'};
            model = model.registerVarNames(varnames);

            if model.use_thermal
                varnames = {'jFace', ...
                            'jFaceBc'};
                model = model.registerVarNames(varnames);                
            end
            
            fn = @ElectronicComponent.updateCurrent;
            inputnames = {'phi', 'conductivity'};
            model = model.registerPropFunction({'j', fn, inputnames});

            fn = @ElectronicComponent.updateChargeConservation;
            inputnames = {'j', 'jBcSource', 'eSource'};
            model = model.registerPropFunction({'chargeCons', fn, inputnames});
            
            fn = @ElectronicComponent.updateConductivity;
            inputnames = {};
            model = model.registerPropFunction({'conductivity', fn, inputnames});
            
            if model.use_thermal
                
                fn = @ElectronicComponent.updateFaceCurrent;
                inputnames = {'j', 'jFaceBc'};
                model = model.registerPropFunction({'jFace', fn, inputnames});

                fn = @ElectronicComponent.updateFaceBcCurrent;
                inputnames = {};
                model = model.registerPropFunction({'jFaceBc', fn, inputnames});
                
            end

        end

        function state = updateConductivity(model, state)
            % default function to update conductivity
            state.conductivity = model.EffectiveElectricalConductivity;
            
        end
        
        
        function state = updateFaceCurrent(model, state)
            
            G = model.G;
            nf = G.faces.num;
            intfaces = model.operators.internalConn;
            
            j       = state.j;
            jFaceBc = state.jFaceBc;
            
            zeroFaceAD = model.AutoDiffBackend.convertToAD(zeros(nf, 1), j);
            jFace = zeroFaceAD + jFaceBc;
            jFace(intfaces) = j;
            
            state.jFace = jFace;
            
        end

        function state = updateFaceBcCurrent(model, state)
            
            state.jFaceBc = 0;
            
        end
        
        function state = updateCurrent(model, state)
        % Assemble electrical current which is stored in :code:`state.j` 

            sigma = state.conductivity;
            phi   = state.phi;
            
            j = assembleFlux(model, phi, sigma); 

            state.j = j;
            
        end
        
        function state = updateChargeConservation(model, state)
        % Assemble residual of the charge conservation equation which is stored in :code:`state.chargeCons`
           
            state = model.updateCurrent(state);

            flux   = state.j;
            bcsource = state.jBcSource;
            source = state.eSource;
            accum  = zeros(model.G.cells.num,1);
            
            chargeCons = assembleConservationEquation(model, flux, bcsource, source, accum);
            
            state.chargeCons = chargeCons;
            
        end
        
    end
end




%{
Copyright 2021-2022 SINTEF Industry, Sustainable Energy Technology
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
