classdef Electrode < BaseModel
%
% The Electrode model is made of two sub-models : an electrode active component (see
% :class:`Electrochemistry.ActiveMaterial`) and a current collector (see
% :class:`Electrochemistry.CurrentCollector`)
%
    properties
        
        ActiveMaterial   % instance of :class:`Electrochemistry.ActiveMaterial`
        CurrentCollector % instance of :class:`Electrochemistry.CurrentCollector`

        couplingTerm
        
        include_current_collector
        
    end

    methods
        
        function model = Electrode(paramobj, params)
        % paramobj is instance of :class:`Electrochemistry.Electrodes.ElectrodeInputParams`
            model = model@BaseModel();
            
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks', false);
            
            fdnames = {'G', ...
                       'couplingTerm'};
            model = dispatchParams(model, paramobj, fdnames);
            
            % Assign the two components
            model.ActiveMaterial = model.setupActiveMaterial(paramobj.ActiveMaterial);
            
            if params.include_current_collector
                assert(~isempty(paramobj.CurrentCollector), 'current collector input data is missing')
                model.CurrentCollector = model.setupCurrentCollector(paramobj.CurrentCollector);
                model.include_current_collector = true;
            else
                assert(isempty(paramobj.CurrentCollector), 'current collector data (and probably grid) is given. Not using it, as required by input flag, will create problem')
                model.include_current_collector = false;
            end
            
        end
        
        function model = registerVarAndPropfuncNames(model)
            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)
            
            model = registerVarAndPropfuncNames@BaseModel(model);
            
            % define shorthands
            am = 'ActiveMaterial';
            cc = 'CurrentCollector';
            am = 'ActiveMaterial';
            
            model = model.registerVarName('T');
            
            fn = @Electrode.updateCoupling;
            inputnames = {{am, 'phi'}, ...
                          {cc , 'phi'}};
            model = model.registerPropFunction({{am, 'jCoupling'}, fn, inputnames});
            model = model.registerPropFunction({{cc , 'jCoupling'}, fn, inputnames});
            model = model.registerPropFunction({{cc , 'eSource'}  , fn, inputnames});
            
            % Temperature coupling between current collector and electrode active component
            inputnames = {{am, 'T'}, ...
                          {cc , 'T'}};
            fn = @Electrode.updateTemperatureCoupling;
            model = model.registerPropFunction({{am, 'jHeatBcSource'}, fn, inputnames});
            model = model.registerPropFunction({{cc , 'jHeatBcSource'}, fn, inputnames});
            
        end
        
        function am = setupActiveMaterial(model, paramobj)
        % paramobj is instance of ActiveMaterialInputParams
        % standard instantiation (ActiveMaterial is specified in ActiveMaterial instantiation)
            am = ActiveMaterial(paramobj);
        end
        
        function cc = setupCurrentCollector(model, paramobj)
        % standard instantiation 
            cc = CurrentCollector(paramobj);
        end
        
        function state = updateCoupling(model, state)
        % setup coupling terms between the current collector and the electrode active component            
            
            if model.include_current_collector
                
                elde  = model;
                am   = 'ActiveMaterial';
                cc    = 'CurrentCollector';

                am_phi = state.(am).phi;
                cc_phi = state.(cc).phi;

                am_sigmaeff = elde.(am).EffectiveElectricalConductivity;
                cc_sigmaeff = elde.(cc).EffectiveElectricalConductivity;
                
                %% We setup the current transfers between CurrentCollector and ActiveMaterial
                
                am_jCoupling  = am_phi*0.0; %NB hack to initialize zero ad
                cc_jCoupling = cc_phi*0.0; %NB hack to initialize zero ad

                coupterm = model.couplingTerm;
                face_cc = coupterm.couplingfaces(:, 1);
                face_am = coupterm.couplingfaces(:, 2);
                [teac, bccell_am] = elde.(am).operators.harmFaceBC(am_sigmaeff, face_am);
                [tcc, bccell_cc] = elde.(cc).operators.harmFaceBC(cc_sigmaeff, face_cc);

                bcphi_am = am_phi(bccell_am);
                bcphi_cc = cc_phi(bccell_cc);

                trans = 1./(1./teac + 1./tcc);
                crosscurrent = trans.*(bcphi_cc - bcphi_am);
                am_jCoupling = subsasgnAD(am_jCoupling,bccell_am,crosscurrent);
                cc_jCoupling = subsasgnAD(cc_jCoupling,bccell_cc, -crosscurrent);

                G = model.(cc).G;
                nf = G.faces.num;
                sgn = model.(cc).operators.sgn;
                zeroFaceAD = model.AutoDiffBackend.convertToAD(zeros(nf, 1), cc_phi);
                cc_jFaceCoupling = zeroFaceAD;
                %cc_jFaceCoupling(face_cc) = sgn(face_cc).*crosscurrent;
                cc_jFaceCoupling = subsasgnAD(cc_jFaceCoupling,face_cc,sgn(face_cc).*crosscurrent);
                assert(~any(isnan(sgn(face_cc))));
                
                G = model.(am).G;
                nf = G.faces.num;
                sgn = model.(am).operators.sgn;
                zeroFaceAD = model.AutoDiffBackend.convertToAD(zeros(nf, 1), am_phi);
                am_jFaceCoupling = zeroFaceAD;
                %am_jFaceCoupling(face_am) = -sgn(face_am).*crosscurrent;
                am_jFaceCoupling = subsasgnAD(am_jFaceCoupling,face_am,-sgn(face_am).*crosscurrent);
                assert(~any(isnan(sgn(face_am))));
                
                % We set here volumetric current source to zero for current collector (could have been done at a more logical place but
                % let us do it here, for simplicity)
                state.(cc).eSource = zeros(elde.(cc).G.cells.num, 1);
                
                state.(am).jCoupling = am_jCoupling;
                state.(cc).jCoupling  = cc_jCoupling;
                
                state.(am).jFaceCoupling = am_jFaceCoupling;
                state.(cc).jFaceCoupling = cc_jFaceCoupling;
                
            end
                
           %     state.(am).jCoupling     = 0;
           %     state.(am).jFaceCoupling = 0;
                
           % end
            
            
        end

        
        function state = updateTemperatureCoupling(model, state)
        % setup coupling terms between the current collector and the electrode active component            
            
            elde  = model;
            am = 'ActiveMaterial';
            cc = 'CurrentCollector';

            am_T = state.(am).T;
            cc_T = state.(cc).T;

            am_tC = elde.(am).thermalConductivity;
            cc_tC = elde.(cc).thermalConductivity;
    
            %% We setup the current transfers between CurrentCollector and ActiveMaterial
            
            am_heatCoupling  = am_T*0.0; %NB hack to initialize zero ad
            cc_heatCoupling = cc_T*0.0; %NB hack to initialize zero ad

            coupterm = model.couplingTerm;
            face_cc = coupterm.couplingfaces(:, 1);
            face_am = coupterm.couplingfaces(:, 2);
            [teac, bccell_am] = elde.(am).operators.harmFaceBC(am_tC, face_am);
            [tcc, bccell_cc] = elde.(cc).operators.harmFaceBC(cc_tC, face_cc);

            bcT_am = am_T(bccell_am);
            bcT_cc = cc_T(bccell_cc);

            trans = 1./(1./teac + 1./tcc);
            crossFluxT = trans.*(bcT_cc - bcT_am);
            am_heatCoupling(bccell_am) = crossFluxT;
            cc_heatCoupling(bccell_cc) = - crossFluxT;

            state.(am).jHeatBcSource = am_heatCoupling;
            state.(cc).jHeatBcSource = cc_heatCoupling;

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
