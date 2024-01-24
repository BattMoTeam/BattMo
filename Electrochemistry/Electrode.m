classdef Electrode < BaseModel
%
% The Electrode model is made of two sub-models : a coating (see
% :class:`Electrochemistry.Coating`) and a current collector (see
% :class:`Electrochemistry.CurrentCollector`)
%
    properties

        %% Sub-Models

        Coating          % instance of :class:`Electrochemistry.ActiveMaterial`
        CurrentCollector % instance of :class:`Electrochemistry.CurrentCollector`

        %% Coupling parameters

        couplingTerm

        %% Computed parameters at setup

        include_current_collectors
        use_thermal

    end

    methods

        function model = Electrode(inputparams)
        % inputparams is instance of :class:`Electrochemistry.Electrodes.ElectrodeInputParams`

            model = model@BaseModel();

            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks', true);

            fdnames = {'G'                         , ...
                       'couplingTerm'              , ...
                       'include_current_collectors', ...
                       'use_thermal'};
            model = dispatchParams(model, inputparams, fdnames);

            model.Coating = Coating(inputparams.Coating);

            if inputparams.include_current_collectors
                model.include_current_collectors = true;
                assert(~isempty(inputparams.CurrentCollector), 'current collector input data is missing')
                model.CurrentCollector = model.setupCurrentCollector(inputparams.CurrentCollector);
            else
                model.include_current_collectors = false;
                % if isempty(inputparams.CurrentCollector.G)
                %    warning('current collector data is given, but we are not using it, as required by input flag');
                % end
            end

        end

        function model = registerVarAndPropfuncNames(model)

        %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)

            % define shorthands
            cc = 'CurrentCollector';
            co = 'Coating';

            if ~model.include_current_collectors
                model.subModelNameList = {co};
            end

            model = registerVarAndPropfuncNames@BaseModel(model);

            if ~model.include_current_collectors
                model = model.registerVarName({co, 'jExternal'});
            end


            if model.include_current_collectors

                fn = @Electrode.updateCoupling;
                inputnames = {{co, 'phi'}, ...
                              {cc , 'phi'}};
                model = model.registerPropFunction({{co, 'jCoupling'}, fn, inputnames});
                model = model.registerPropFunction({{cc, 'jCoupling'}, fn, inputnames});
                model = model.registerPropFunction({{cc, 'eSource'}  , fn, inputnames});

                if model.use_thermal
                    model = model.registerPropFunction({{co, 'jFaceCoupling'}, fn, inputnames});
                    model = model.registerPropFunction({{cc, 'jFaceCoupling'}, fn, inputnames});
                end

            end

        end

        function cc = setupCurrentCollector(model, inputparams)
        % standard instantiation
            cc = CurrentCollector(inputparams);
        end

        function model = setTPFVgeometry(model, tPFVgeometry)
        % tPFVgeometry should be instance of TwoPointFiniteVolumeGeometry or MutableTwoPointFiniteVolumeGeometry

            co = 'Coating';

            model.G.parentGrid.tPFVgeometry = tPFVgeometry;

            model.(co) = model.(co).setTPFVgeometry(tPFVgeometry);

            if model.include_current_collectors
                cc = 'CurrentCollector';
                model.(cc) = model.(cc).setTPFVgeometry(tPFVgeometry);
            end

        end

        function state = updateCoupling(model, state)
        % setup coupling terms between the current collector and the electrode active component

            if model.include_current_collectors

                elde  = model;

                co = 'Coating';
                cc = 'CurrentCollector';

                co_phi = state.(co).phi;
                cc_phi = state.(cc).phi;

                co_sigmaeff = elde.(co).effectiveElectronicConductivity;
                cc_sigmaeff = elde.(cc).effectiveElectronicConductivity;

                %% We setup the current transfers between CurrentCollector and ActiveMaterial

                co_jCoupling  = co_phi*0.0; %NB hack to initialize zero ad
                cc_jCoupling = cc_phi*0.0; %NB hack to initialize zero ad

                coupterm = model.couplingTerm;
                face_cc = coupterm.couplingfaces(:, 1);
                face_co = coupterm.couplingfaces(:, 2);
                [teac, bccell_co] = elde.(co).operators.transFaceBC(face_co);
                teac = co_sigmaeff*teac;
                [tcc, bccell_cc] = elde.(cc).operators.transFaceBC( face_cc);
                tcc = cc_sigmaeff*teac;

                bcphi_co = co_phi(bccell_co);
                bcphi_cc = cc_phi(bccell_cc);

                trans = 1./(1./teac + 1./tcc); % Harmonic average
                crosscurrent = trans.*(bcphi_cc - bcphi_co);
                co_jCoupling = subsasgnAD(co_jCoupling,bccell_co, crosscurrent);
                cc_jCoupling = subsasgnAD(cc_jCoupling,bccell_cc, -crosscurrent);

                G = model.(cc).G;
                nf = G.getNumberOfFaces();
                sgn = model.(cc).operators.sgn;
                zeroFaceAD = model.AutoDiffBackend.convertToAD(zeros(nf, 1), cc_phi);
                cc_jFaceCoupling = zeroFaceAD;
                cc_jFaceCoupling = subsasgnAD(cc_jFaceCoupling,face_cc, sgn(face_cc).*crosscurrent);
                assert(~any(isnan(sgn(face_cc))));

                G = model.(co).G;
                nf = G.getNumberOfFaces();
                sgn = model.(co).operators.sgn;
                zeroFaceAD = model.AutoDiffBackend.convertToAD(zeros(nf, 1), co_phi);
                co_jFaceCoupling = zeroFaceAD;
                co_jFaceCoupling = subsasgnAD(co_jFaceCoupling, face_co, -sgn(face_co).*crosscurrent);
                assert(~any(isnan(sgn(face_co))));

                % We set here volumetric current source to zero for current collector (could have been done at a more
                % logical place but let us do it here, for simplicity)
                state.(cc).eSource = zeros(elde.(cc).G.getNumberOfCells(), 1);

                state.(co).jCoupling = co_jCoupling;
                state.(cc).jCoupling = cc_jCoupling;

                state.(co).jFaceCoupling = co_jFaceCoupling;
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
