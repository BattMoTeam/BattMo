classdef KineticParamSetter
%%
% For this parameter setter class, the chosen parameters are one for each rate limited (kinetic) process
%
%  - Chemical reaction                        : Volumetric surface areas
%  - Charge transport in electrodes           : Bruggeman coefficient
%  - Diffusion in solid electrodes            : Diffusion coefficient
%  - Mass and charge transport in electrolyte : Bruggeman coefficient
%
    
    properties

        boxLims

    end

    methods

        function paramsetter = KineticParamSetter()

            % Some default values for the bounding box are given here but they should be checked.
            boxLims = [[1e4, 1e9];
                       [1e4, 1e9];
                       [0.5, 3];
                       [0.5, 3];
                       [1e-14, 1e-9];
                       [1e-14, 1e-9];
                       [0.5, 3];
                      ];

            paramsetter.boxLims = boxLims;

        end

        function simsetup = setValues(paramsetter, simsetup, X)
        %%
        % Given the vector of parameters X, set the model with those parameters

            model = simsetup.model;
            
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            elyte   = 'Electrolyte';
            co      = 'Coating';
            am      = 'ActiveMaterial';
            sd      = 'SolidDiffusion';
            itf     = 'Interface';

            %% Automatic part of the update
            locs = paramsetter.locations();

            for k = 1:numel(locs)
                loc = locs{k};
                model = setfield(model, loc{:}, X(k));
            end

            %% Update dependencies (manual work)
            eldes = {ne, pe};

            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};

                bg    = model.(elde).(co).bruggemanCoefficient;
                kappa = model.(elde).(co).electronicConductivity;
                vf    = model.(elde).(co).volumeFraction;

                model.(elde).(co).effectiveElectronicConductivity = kappa*vf^bg;

                model.(elde).(co).(am).(sd).volumetricSurfaceArea = model.(elde).(co).(am).(itf).volumetricSurfaceArea;

            end

            if model.(elyte).useRegionBruggemanCoefficients

                nc = model.(elyte).G.getNumberOfCells();
                bg = zeros(nc, 1);

                tagmap = struct('NegativeElectrode', 1, 'PositiveElectrode', 2, 'Separator', 3);
                tags = model.(elyte).regionTags;

                for k = 1:numel(locs)
                    loc = locs{k};
                    if strcmp(loc{1}, elyte) && strcmp(loc{2}, 'regionBruggemanCoefficients')
                        % Should be set by the automatic update above:
                        % bval = X(k);
                        % model.(elyte).regionBruggemanCoefficients.(loc{3}) = bval;
                        bval = model.(elyte).regionBruggemanCoefficients.(loc{3});
                        bg = subsetPlus(bg, bval, (tags == tagmap.(loc{3})));
                    end

                end

                model.(elyte).bruggemanCoefficient = bg;

            end

            model.jsonstruct = [];

            simsetup.model = model;
            
        end


        function X = getValues(paramsetter, simsetup)
        %%
        % Given a model, retrieve the value of the optimization parameter in a vector X

            model = simsetup.model;
            
            locs  = paramsetter.locations();
            short = paramsetter.shortlocs();

            vals = struct();

            for k = 1:numel(locs)
                loc = locs{k};
                mval = getfield(model, loc{:});
                s = short{k};
                vals = setfield(vals, s{:}, mval);
            end

            X = paramsetter.setToVector(vals);

        end

        function params = setupModelParameters(paramsetter, simsetup, varargin)

            function vals = localGetValues(model)

                newsimsetup = simsetup;
                newsimsetup.model = model;
                vals = paramsetter.getValues(newsimsetup);
                
            end

            function model = localSetValues(model, v)

                newsimsetup = simsetup;
                newsimsetup.model = model;
                newsimsetup = paramsetter.setValues(newsimsetup, v);
                model = newsimsetup.model;
                
            end

            getValues = @(model, notused) localGetValues(model);
            setValues = @(model, notused, v) localSetValues(model, v);
            
            params{1} = ModelParameter(simsetup, ...
                                       'name'     , 'ParamSetter'      , ...
                                       'belongsTo', 'model'            , ...
                                       'location' , {''}               , ...
                                       'boxLims'  , paramsetter.boxLims, ...
                                       'getfun'   , getValues          , ...
                                       'setfun'   , setValues          , ...
                                       varargin{:});
            

        end
        
        function vals = setFromVector(paramsetter, X)
        %%
        % Conversion from vector format to structure with human readable names
            slocs = paramsetter.shortlocs();
            vals = struct();

            for k = 1:numel(slocs)
                s = slocs{k};
                vals = setfield(vals, s{:}, X(k));
            end

        end

        function X = setToVector(paramsetter, vals)
        %%
        % Conversion from structure with human readable names to vector format
            locs = paramsetter.shortlocs();
            X = nan(numel(locs), 1);

            for k = 1:numel(locs)
                loc = locs{k};
                X(k) = getfield(vals, loc{:});
            end

        end

        function locs = locations(paramsetter)
        % Utility function, which gives the name of where the values of the parameter are stored in the model
            
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            co    = 'Coating';
            cc    = 'CurrentCollector';
            am    = 'ActiveMaterial';
            sd    = 'SolidDiffusion';
            itf   = 'Interface';
            elyte = 'Electrolyte';
            sep   = 'Separator';

            locs = {};

            locs{1} = {ne, co, am, itf, 'volumetricSurfaceArea'};
            locs{2} = {pe, co, am, itf, 'volumetricSurfaceArea'};

            locs{3} = {ne, co, 'bruggemanCoefficient'};
            locs{4} = {pe, co, 'bruggemanCoefficient'};

            locs{5} = {ne, co, am, sd, 'referenceDiffusionCoefficient'};
            locs{6} = {pe, co, am, sd, 'referenceDiffusionCoefficient'};

            locs{7} = {elyte, 'bruggemanCoefficient'};

            assert(numel(locs) == size(paramsetter.boxLims,1), 'Number of locations and box limits must match');

        end

        function s = shortlocs(paramsetter)
        % Utility function that provides short names for the optimization parameters
            s = cellfun(@(loc) {loc{1}, loc{end}}, paramsetter.locations(), 'un', false);

        end

        function sn = shortnames(paramsetter, join)
        % Utility function that provides the conversion from long to short name
            s = cellfun(@(loc) {loc{1}, loc{end}}, paramsetter.locations(), 'un', false);

            if nargin < 2
                join = false;
            end

            dict = struct('NegativeElectrode'              , 'ne'   , ...
                          'PositiveElectrode'              , 'pe'   , ...
                          'effectiveElectronicConductivity', 'kappa', ...
                          'Electrolyte'                    , 'elyte', ...
                          'volumetricSurfaceArea'          , 'vsa'  , ...
                          'bruggemanCoefficient'           , 'bg'   , ...
                          'referenceDiffusionCoefficient'  , 'D'    , ...
                          'Separator'                      , 'sep');

            slocs = paramsetter.shortlocs();
            sn = cell(numel(slocs), 1);

            for k = 1:numel(slocs)
                loc = slocs{k};
                sn{k} = [dict.(loc{1}), '_', dict.(loc{2})];
            end

            if join
                sn = strjoin(sn, ' ');
            end

        end

        function printBoxLims(paramsetter)
            
            boxLims = paramsetter.boxLims;
            boxLims = num2cell(boxLims);

            params = paramsetter.shortnames();

            T = cell2table(horzcat(params, boxLims), 'VariableNames', {'parameters', 'min', 'max'});
            display(T);
            
        end
        
        function paramsetter = setupRelBox(paramsetter, relsize, X)

            if numel(relsize) == 1
                relsize = relsize*ones(size(X));
            end

            boxLims = paramsetter.boxLims;
            for k = 1:numel(X)
                boxLims(k, 1) = X(k)*(1 - relsize(k));
                boxLims(k, 2) = X(k)*(1 + relsize(k));
            end

            paramsetter.boxLims = boxLims;

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
