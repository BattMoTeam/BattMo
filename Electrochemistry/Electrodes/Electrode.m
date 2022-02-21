classdef Electrode < BaseModel
%
% The Electrode model is made of two sub-models : an electrode active component (see
% :class:`Electrochemistry.ElectrodeActiveComponent`) and a current collector (see
% :class:`Electrochemistry.CurrentCollector`)
%
    properties
        
        ElectrodeActiveComponent % instance of :class:`Electrochemistry.ElectrodeActiveComponent`
        CurrentCollector         % instance of :class:`Electrochemistry.CurrentCollector`
        couplingTerm
        
    end

    methods
        
        function model = Electrode(paramobj)
        % paramobj is instance of :class:`Electrochemistry.Electrodes.ElectrodeInputParams`
            model = model@BaseModel();
            
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks', false);
            
            fdnames = {'G', ...
                       'couplingTerm'};
            model = dispatchParams(model, paramobj, fdnames);
            
            % Assign the two components
            model.ElectrodeActiveComponent = model.setupElectrodeActiveComponent(paramobj.ElectrodeActiveComponent);
            model.CurrentCollector = model.setupCurrentCollector(paramobj.CurrentCollector);

            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)
            
            model = model.registerSubModels({'ElectrodeActiveComponent', 'CurrentCollector'});
            
            model = model.registerVarName('T');
            
            fn = @Electrode.updateCoupling;
            inputnames = {VarName({'eac'}, 'phi'), ...
                          VarName({'cc'}, 'phi')};
            model = model.registerPropFunction({VarName({'eac'}, 'jCoupling'), fn, inputnames});
            model = model.registerPropFunction({VarName({'cc'} , 'jCoupling'), fn, inputnames});
            model = model.registerPropFunction({VarName({'cc'} , 'eSource')  , fn, inputnames});
            
            % Temperature coupling between current collector and electrode active component
            inputnames = {VarName({'eac'}, 'T'), ...
                          VarName({'cc'}, 'T')};
            fn = @Electrode.updateTemperatureCoupling;
            model = model.registerPropFunction({VarName({'eac'}, 'jHeatBcSource'), fn, inputnames});
            model = model.registerPropFunction({VarName({'cc'}, 'jHeatBcSource'), fn, inputnames});
            
        end
        
        function eac = setupElectrodeActiveComponent(model, paramobj)
        % paramobj is instanceo of ElectrodeActiveComponentInputParams
        % standard instantiation (ActiveMaterial is specified in ElectrodeActiveComponent instantiation)
            eac = ElectrodeActiveComponent(paramobj);
        end
        
        function cc = setupCurrentCollector(model, paramobj)
        % standard instantiation 
            cc = CurrentCollector(paramobj);
        end
        
        function state = updateCoupling(model, state)
        % setup coupling terms between the current collector and the electrode active component            
            
            elde  = model;
            eac   = 'ElectrodeActiveComponent';
            cc    = 'CurrentCollector';

            eac_phi = state.(eac).phi;
            cc_phi = state.(cc).phi;

            eac_sigmaeff = elde.(eac).EffectiveElectricalConductivity;
            cc_sigmaeff = elde.(cc).EffectiveElectricalConductivity;
    
            %% We setup the current transfers between CurrentCollector and ElectrodeActiveComponent
            
            eac_jCoupling  = eac_phi*0.0; %NB hack to initialize zero ad
            cc_jCoupling = cc_phi*0.0; %NB hack to initialize zero ad

            coupterm = model.couplingTerm;
            face_cc = coupterm.couplingfaces(:, 1);
            face_eac = coupterm.couplingfaces(:, 2);
            [teac, bccell_eac] = elde.(eac).operators.harmFaceBC(eac_sigmaeff, face_eac);
            [tcc, bccell_cc] = elde.(cc).operators.harmFaceBC(cc_sigmaeff, face_cc);

            bcphi_eac = eac_phi(bccell_eac);
            bcphi_cc = cc_phi(bccell_cc);

            trans = 1./(1./teac + 1./tcc);
            crosscurrent = trans.*(bcphi_cc - bcphi_eac);
            eac_jCoupling(bccell_eac) = crosscurrent;
            cc_jCoupling(bccell_cc) = -crosscurrent;

            G = model.(cc).G;
            nf = G.faces.num;
            sgn = model.(cc).operators.sgn;
            zeroFaceAD = model.AutoDiffBackend.convertToAD(zeros(nf, 1), cc_phi);
            cc_jFaceCoupling = zeroFaceAD;
            cc_jFaceCoupling(face_cc) = sgn(face_cc).*crosscurrent;
            assert(~any(isnan(sgn(face_cc))));
            
            G = model.(eac).G;
            nf = G.faces.num;
            sgn = model.(eac).operators.sgn;
            zeroFaceAD = model.AutoDiffBackend.convertToAD(zeros(nf, 1), eac_phi);
            eac_jFaceCoupling = zeroFaceAD;
            eac_jFaceCoupling(face_eac) = -sgn(face_eac).*crosscurrent;
            assert(~any(isnan(sgn(face_eac))));
            
            % We set here volumetric current source to zero for current collector (could have been done at a more logical place but
            % let us do it here, for simplicity)
            state.(cc).eSource = zeros(elde.(cc).G.cells.num, 1);
            
            state.(eac).jCoupling = eac_jCoupling;
            state.(cc).jCoupling  = cc_jCoupling;
            
            state.(eac).jFaceCoupling = eac_jFaceCoupling;
            state.(cc).jFaceCoupling  = cc_jFaceCoupling;
            
        end

        
        function state = updateTemperatureCoupling(model, state)
        % setup coupling terms between the current collector and the electrode active component            
            
            elde  = model;
            eac   = 'ElectrodeActiveComponent';
            cc    = 'CurrentCollector';

            eac_T = state.(eac).T;
            cc_T = state.(cc).T;

            eac_tC = elde.(eac).thermalConductivity;
            cc_tC = elde.(cc).thermalConductivity;
    
            %% We setup the current transfers between CurrentCollector and ElectrodeActiveComponent
            
            eac_heatCoupling  = eac_T*0.0; %NB hack to initialize zero ad
            cc_heatCoupling = cc_T*0.0; %NB hack to initialize zero ad

            coupterm = model.couplingTerm;
            face_cc = coupterm.couplingfaces(:, 1);
            face_eac = coupterm.couplingfaces(:, 2);
            [teac, bccell_eac] = elde.(eac).operators.harmFaceBC(eac_tC, face_eac);
            [tcc, bccell_cc] = elde.(cc).operators.harmFaceBC(cc_tC, face_cc);

            bcT_eac = eac_T(bccell_eac);
            bcT_cc = cc_T(bccell_cc);

            trans = 1./(1./teac + 1./tcc);
            crossFluxT = trans.*(bcT_cc - bcT_eac);
            eac_heatCoupling(bccell_eac) = crossFluxT;
            cc_heatCoupling(bccell_cc) = - crossFluxT;

            state.(eac).jHeatBcSource = eac_heatCoupling;
            state.(cc).jHeatBcSource = cc_heatCoupling;

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
