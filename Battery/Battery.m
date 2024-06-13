classdef Battery < GenericBattery
%
% The battery model consists of
%
% * an Electrolyte model given in :attr:`Electrolyte` property
% * a Negative Electrode Model given in :attr:`NegativeElectrode` property
% * a Positive Electrode Model given in :attr:`PositiveElectrode` property
% * a Thermal model given in :attr:`ThermalModel` property
%
    properties

        equationIndices

    end

    methods

        function model = Battery(inputparams)

            model = model@GenericBattery(inputparams);

            % setup equations and variable names selected in the model
            model = model.setupSelectedModel();

        end

        function model = setupSelectedModel(model, varargin)
        % The system of equation should fullfill a special structure to fit into the iterative linear solver with
        % preconditioner. We create this structure here.

            opt = struct('reduction', []);
            opt = merge_options(opt, varargin{:});
            % For the reduction structure format, see battmodDir()/Utilities/JsonSchemas/linearsolver.schema.json and the reduction property

            elyte   = 'Electrolyte';
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            co      = 'Coating';
            am      = 'ActiveMaterial';
            am1     = 'ActiveMaterial1';
            am2     = 'ActiveMaterial2';
            sd      = 'SolidDiffusion';
            cc      = 'CurrentCollector';
            ctrl    = 'Control';
            thermal = 'ThermalModel';

            addedVarNames = {};

            varEqTypes ={{elyte, 'c'}   , {elyte,  'massCons'}     , 'cell'; ...
                         {elyte, 'phi'} , {elyte,  'chargeCons'}   , 'cell'; ...
                         {ne, co, 'phi'}, {ne, co, 'chargeCons'}   , 'cell'; ...
                         {pe, co, 'phi'}, {pe, co, 'chargeCons'}   , 'cell'; ...
                         {ctrl, 'E'}    , {ctrl, 'EIequation'}     , 'ctrl'; ...
                         {ctrl, 'I'}    , {ctrl, 'controlEquation'}, 'ctrl'};

            eldes = {ne, pe};


            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};

                if model.(elde).(co).activeMaterialModelSetup.composite
                    ams = {am1, am2};
                else
                    ams = {am};
                end

                for iam = 1 : numel(ams)
                    amc = ams{iam};
                    switch model.(elde).(co).(amc).diffusionModelType
                      case 'simple'
                        newentries = {{elde, co, amc, sd, 'cAverage'}, {elde, co, amc, sd, 'massCons'}        , 'cell'; ...
                                      {elde, co, amc, sd, 'cSurface'}, {elde, co, amc, sd, 'solidDiffusionEq'}, 'cell'};
                      case 'full'
                        newentries = {{elde, co, amc, sd, 'c'}       , {elde, co, amc, sd, 'massCons'}        , 'cell'; ...
                                      {elde, co, amc, sd, 'cSurface'}, {elde, co, amc, sd, 'solidDiffusionEq'}, 'cell'};
                      otherwise
                        error('diffusionModelType not recognized');
                    end

                    varEqTypes = vertcat(varEqTypes, newentries);

                end

                if model.include_current_collectors

                    newentries = {{elde, cc, 'phi'}, {elde, cc, 'chargeCons'}, 'cell'};
                    varEqTypes = vertcat(varEqTypes, newentries);

                end

            end

            if model.use_thermal

                newentries = {{thermal, 'T'}, {thermal, 'energyCons'}, 'cell'};
                varEqTypes = vertcat(varEqTypes, newentries);

            end

            primaryVarNames = varEqTypes(:, 1);
            equationTypes   = varEqTypes(:, 3);

            % The variable and equation lists are not a priori ordered (in the sense that we have 'cell' types first and
            % the other types after). It is a requirement in some setup of the linear solver.
            % Note : if you use a direct solver, this is not used.

            variableReordered = false;

            if ~isempty(opt.reduction)
                reduc = opt.reduction;
                % We set the type of the variable to be reduced as 'reduced' or 'specialReduced' (anything but 'cell')
                % and we move them at the end of list respecting the order in which they have been given.
                if ~isempty(reduc) && reduc.doReduction
                    neq = numel(equationTypes);
                    equationTypes = cell(neq, 1);
                    for ieqtype = 1 : neq
                        equationTypes{ieqtype} = 'cell';
                    end

                    variables = reduc.variables;
                    if isstruct(variables)
                        variables = num2cell(variables);
                    end
                    einds = nan(numel(variables),1);
                    for ivar = 1 : numel(variables)
                        var = variables{ivar};
                        [found, ind] = Battery.getVarIndex(var.name, primaryVarNames);
                        if ~found
                            error('variable to be reduce has not been found');
                        end
                        equationTypes{ind} = 'reduced';
                        if isfield(var, "special") && var.special
                            equationTypes{ind} = 'specialReduced';
                        end
                        einds(ivar) = ind;
                        order(ivar) = var.order;
                    end

                    [~, ind] = sort(order);
                    einds = einds(ind);

                    inds = (1 : neq)';
                    inds(einds) = [];
                    inds = [inds; einds];
                    variableReordered = true;

                end
            end

            if ~variableReordered
                % We reorder to get first the 'cells' type (required for reduction in iterative solver)
                iscell = ismember(equationTypes, {'cell'});
                inds = [find(iscell); find(~iscell)];
            end


            primaryVarNames  = varEqTypes(inds, 1);
            equationVarNames = varEqTypes(inds, 2);
            equationTypes    = equationTypes(inds);

            % We use shortened names for easier visualisation and because Matlab also has a limitation on the lenght of
            % a field name in a structure.
            function str = setupName(varname)
                shortvarname = cellfun(@(elt) Battery.shortenName(elt), varname, 'uniformoutput', false);
                str = Battery.varToStr(shortvarname);
            end
            equationNames = cellfun(@(varname) setupName(varname), equationVarNames, 'uniformoutput', false);

            equationIndices = struct();
            for ieq = 1 : numel(equationNames)
                equationIndices.(equationNames{ieq}) = ieq;
            end

            model.primaryVarNames  = primaryVarNames;
            model.equationVarNames = equationVarNames;
            model.equationNames    = equationNames;
            model.equationTypes    = equationTypes;
            model.equationIndices  = equationIndices;

        end


        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
        % Assembly of the governing equation

            opts = struct('ResOnly', false, 'iteration', 0, 'reverseMode', false);
            opts = merge_options(opts, varargin{:});

            time = state0.time + dt;

            if (not(opts.ResOnly) && not(opts.reverseMode))
                state = model.initStateAD(state);
            elseif opts.reverseMode
               dispif(mrstVerbose, 'No AD initialization in equation old style')
               state0 = model.initStateAD(state0);
            else
                assert(opts.ResOnly);
            end


            %% We call the assembly equations ordered from the graph

            funcCallList = model.funcCallList;

            for ifunc = 1 : numel(funcCallList)
                eval(funcCallList{ifunc});
            end

            %% We apply some scaling

            % Shorthands used in this function
            battery = model;
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            co      = 'Coating';
            am      = 'ActiveMaterial';
            am1     = 'ActiveMaterial1';
            am2     = 'ActiveMaterial2';
            elyte   = 'Electrolyte';
            am      = 'ActiveMaterial';
            itf     = 'Interface';
            sd      = 'SolidDiffusion';
            thermal = 'ThermalModel';
            ctrl    = 'Control';

            eldes = {ne, pe};

            massConsScaling = model.con.F;

            state.(elyte).massCons = state.(elyte).massCons.*massConsScaling;

            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};

                if model.(elde).(co).activeMaterialModelSetup.composite
                    ams = {am1, am2};
                else
                    ams = {am};
                end

                for iam = 1 : numel(ams)

                    amc = ams{iam};

                    switch model.(elde).(co).(amc).diffusionModelType

                      case 'simple'

                        state.(elde).(co).(amc).(sd).massCons         = massConsScaling.*state.(elde).(co).(amc).(sd).massCons;
                        state.(elde).(co).(amc).(sd).solidDiffusionEq = massConsScaling.*battery.(elde).(co).G.getVolumes()/dt.*state.(elde).(co).(amc).(sd).solidDiffusionEq;

                      case 'full'

                        n    = model.(elde).(co).(amc).(itf).numberOfElectronsTransferred;
                        F    = model.con.F;
                        vol  = model.(elde).(co).G.getVolumes();
                        rp   = model.(elde).(co).(amc).(sd).particleRadius;
                        vsf  = model.(elde).(co).(amc).(sd).volumetricSurfaceArea;

                        surfp = 4*pi*rp^2;

                        scalingcoef = (vsf*vol(1)*n*F)/surfp;

                        state.(elde).(co).(amc).(sd).massCons         = scalingcoef.*state.(elde).(co).(amc).(sd).massCons;
                        state.(elde).(co).(amc).(sd).solidDiffusionEq = scalingcoef.*state.(elde).(co).(amc).(sd).solidDiffusionEq;

                      otherwise

                        error('diffusionModelType not recognized');

                    end
                end
            end

            for ieq = 1 : numel(model.equationVarNames)
                eqs{ieq} = model.getProp(state, model.equationVarNames{ieq});
            end

            ei = model.equationIndices;

            % By doing this linear transformation, we remove the direct dependency of the mass conservation equation with
            % respect to the potential gradien. Only the concentration gradient remains in the equation. This is a special
            % property when the transference t is a constant.

            mieq = ei.(Battery.varToStr({'elyte', 'massCons'}));
            cieq = ei.(Battery.varToStr({'elyte', 'chargeCons'}));

            t = model.Electrolyte.species.transferenceNumber;

            eqs{mieq} = eqs{mieq} - t*eqs{cieq};

            names       = model.equationNames;
            types       = model.equationTypes;
            primaryVars = model.primaryVarNames;

            %% The equations are reordered in a way that is consitent with the linear iterative solver
            % (the order of the equation does not matter if we do not use an iterative solver)
            ctrltype = state.Control.ctrlType;
            switch ctrltype

              case {'constantCurrent', 'CC_discharge1', 'CC_discharge2', 'CC_charge1', 'charge', 'discharge'}

                eqname = Battery.varToStr({'ctrl', 'EIequation'});
                types{ei.(eqname)} = 'cell';

              case {'constantVoltage', 'CV_charge2'}

                eieqname = Battery.varToStr({'ctrl', 'EIequation'});
                cteqname = Battery.varToStr({'ctrl', 'controlEquation'});

                eqs([ei.(eieqname), ei.(cteqname)]) = eqs([ei.(cteqname), ei.(eieqname)]);

              otherwise

                error('control type not recognized');

            end


            %% Setup LinearizedProblem that can be processed by MRST Newton API
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

        end


    end


end



%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
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
