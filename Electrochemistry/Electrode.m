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

        electrode_case 
        
        include_current_collectors

        use_thermal
        
    end

    methods
        
        function model = Electrode(paramobj)
        % paramobj is instance of :class:`Electrochemistry.Electrodes.ElectrodeInputParams`
            
            model = model@BaseModel();
            
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks', false);
            
            fdnames = {'G', ...
                       'couplingTerm', ...
                       'electrode_case', ...
                       'use_thermal'};
            model = dispatchParams(model, paramobj, fdnames);
            
            switch paramobj.electrode_case
              case 'default'
                model.ActiveMaterial = ActiveMaterial(paramobj.ActiveMaterial);
              case 'composite'
                model.ActiveMaterial = CompositeActiveMaterial(paramobj.ActiveMaterial);
              otherwise
                error('electrode_case not recognized');
            end
            
            if paramobj.include_current_collectors
                model.include_current_collectors = true;
                assert(~isempty(paramobj.CurrentCollector), 'current collector input data is missing')
                model.CurrentCollector = model.setupCurrentCollector(paramobj.CurrentCollector);
            else
                model.include_current_collectors = false;
                % if isempty(paramobj.CurrentCollector.G)
                %    warning('current collector data is given, but we are not using it, as required by input flag');
                % end
            end
                   
        end
        
        function model = registerVarAndPropfuncNames(model)

        %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)

            % define shorthands
            cc = 'CurrentCollector';
            am = 'ActiveMaterial';

            if ~model.include_current_collectors
                model.subModelNameList = {am};
            end
           
            model = registerVarAndPropfuncNames@BaseModel(model);

            if ~model.include_current_collectors
                model = model.registerVarName(VarName({am}, 'jExternal'));
            end
            
            
            if model.include_current_collectors

                fn = @Electrode.updateCoupling;
                inputnames = {{am, 'phi'}, ...
                              {cc , 'phi'}};
                model = model.registerPropFunction({{am, 'jCoupling'}, fn, inputnames});
                model = model.registerPropFunction({{cc, 'jCoupling'}, fn, inputnames});
                model = model.registerPropFunction({{cc, 'eSource'}  , fn, inputnames});

                if model.use_thermal
                    model = model.registerPropFunction({{am, 'jFaceCoupling'}, fn, inputnames});
                    model = model.registerPropFunction({{cc, 'jFaceCoupling'}, fn, inputnames});
                end
                
            end
            
        end
        
        function cc = setupCurrentCollector(model, paramobj)
        % standard instantiation 
            cc = CurrentCollector(paramobj);
        end

        function state = updateCoupling(model, state)
        % setup coupling terms between the current collector and the electrode active component            
            
            if model.include_current_collectors
                
                elde  = model;
                
                am = 'ActiveMaterial';
                cc = 'CurrentCollector';

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
                [teac, bccell_am, bcsgn_am] = elde.(am).G.getBcHarmFace(am_sigmaeff, face_am);
                [tcc, bccell_cc, bcsgn_cc]  = elde.(cc).G.getBcHarmFace(cc_sigmaeff, face_cc);

                bcphi_am = am_phi(bccell_am);
                bcphi_cc = cc_phi(bccell_cc);

                trans = 1./(1./teac + 1./tcc);
                crosscurrent = trans.*(bcphi_cc - bcphi_am);
                am_jCoupling = subsasgnAD(am_jCoupling,bccell_am, crosscurrent);
                cc_jCoupling = subsasgnAD(cc_jCoupling,bccell_cc, -crosscurrent);

                G = model.(cc).G;
                nf = G.topology.faces.num;
                zeroFaceAD = model.AutoDiffBackend.convertToAD(zeros(nf, 1), cc_phi);
                cc_jFaceCoupling = zeroFaceAD;
                cc_jFaceCoupling = subsasgnAD(cc_jFaceCoupling, face_cc, bcsgn_cc.*crosscurrent);
                
                G = model.(am).G; 
                nf = G.topology.faces.num; 
                zeroFaceAD = model.AutoDiffBackend.convertToAD(zeros(nf, 1), am_phi); 
                am_jFaceCoupling = zeroFaceAD; 
                am_jFaceCoupling = subsasgnAD(am_jFaceCoupling, face_am, -bcsgn_am.*crosscurrent); 
                
                % We set here volumetric current source to zero for current collector (could have been done at a more logical place but
                % let us do it here, for simplicity)
                state.(cc).eSource = zeros(elde.(cc).G.getNumberOfCells(), 1);
                
                state.(am).jCoupling = am_jCoupling;
                state.(cc).jCoupling = cc_jCoupling;
                
                state.(am).jFaceCoupling = am_jFaceCoupling;
                state.(cc).jFaceCoupling = cc_jFaceCoupling;
                
            end
                
        end

    end    
end




%{
Copyright 2021-2023 SINTEF Industry, Sustainable Energy Technology
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
