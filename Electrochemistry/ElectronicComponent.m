classdef ElectronicComponent < BaseModel
%
% The ElectronicComponent class is model to assemble charge conservation equation
% 
    properties
        
        EffectiveElectricalConductivity % Effective electrical conductivity
        constants
        
    end

    methods
        
        function model = ElectronicComponent(paramobj)
        % Here, :code:`paramobj` is instance of :class:`ElectronicComponentInputParams <Electrochemistry.ElectronicComponentInputParams>`

            model = model@BaseModel();
            
            % OBS : All the models should have same backend (this is not assigned automaticallly for the moment)
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks', false);
            
            fdnames = {'G', ...
                       'EffectiveElectricalConductivity'};
            
            model = dispatchParams(model, paramobj, fdnames);
            
            % setup discrete differential operators
            model.operators = localSetupOperators(model.G, 'assembleCellFluxOperator', true);
            
            model.constants = PhysicalConstants();

        end

        function model = registerVarAndPropfuncNames(model)
            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)
            
            varnames = {'T'        , ...
                        'phi'      , ...
                        'jBcSource', ...
                        'eSource'  , ...
                        'j'        , ...
                        'chargeCons'};
            model = model.registerVarNames(varnames);


            fn = @ElectronicComponent.updateCurrent;
            inputnames = {'phi'};
            model = model.registerPropFunction({'j', fn, inputnames});
            
            fn = @ElectronicComponent.updateChargeConservation;
            inputnames = {'j', 'jBcSource', 'eSource'};
            model = model.registerPropFunction({'chargeCons', fn, inputnames});
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
        
        
        function state = updateCurrent(model, state)
        % Assemble electrical current which is stored in :code:`state.j` 
            sigmaeff = model.EffectiveElectricalConductivity;
            phi = state.phi;
            
            j = assembleFlux(model, phi, sigmaeff); 

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
Copyright 2009-2021 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BatMo

BatMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BatMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BatMo.  If not, see <http://www.gnu.org/licenses/>.
%}
