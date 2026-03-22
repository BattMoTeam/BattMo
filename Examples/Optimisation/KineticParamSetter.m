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


        active_parameters_inds
        active_parameters_shortnames
    end

    methods

        function paramsetter = KineticParamSetter(varargin)

            opt = struct('active_parameters_inds'                   , []   , ...
                         'active_parameters_shortnames'             , []);
            opt = merge_options(opt, varargin{:});

            
            % Some default values for the bounding box are given here but they should be checked.
            boxLims = [[1e4, 1e7];
                       [1e4, 1e7];
                       [0.5, 3];
                       [0.5, 3];
                       [1e-14, 1e-11];
                       [1e-14, 1e-11];
                       [0.5, 3];
                      ];

            paramsetter.boxLims = boxLims;

            sn = paramsetter.shortnames(paramsetter.allLocations);

            if isempty(opt.active_parameters_inds)

                if isempty(opt.active_parameters_shortnames)
                    active_parameters_inds = (1 : size(boxLims, 1));
                else
                    psn = opt.active_parameters_shortnames;
                    active_parameters_inds = find(ismember(sn, psn));
                end
            else
                active_parameters_inds = opt.active_parameters_inds
            end
            
            paramsetter.boxLims                      = boxLims(active_parameters_inds, :);
            paramsetter.active_parameters_inds       = active_parameters_inds;
            paramsetter.active_parameters_shortnames = sn(active_parameters_inds);
            
        end

        function locs = allLocations(paramsetter)

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

        end

        function locs = locations(paramsetter)

            locs = paramsetter.allLocations();
            locs = locs(paramsetter.active_parameters_inds);

            assert(numel(locs) == size(paramsetter.boxLims, 1), 'Number of locations and box limits must match');

        end

        function s = shortlocs(paramsetter, locations)

            s = cellfun(@(loc) {loc{1}, loc{end}}, locations, 'un', false);

        end

        function sn = shortnames(paramsetter, locations, join)

            if (nargin < 2) | (nargin >= 2 && isempty(locations))
                locations = paramsetter.locations();
            end
            
            if nargin < 3
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

            slocs = paramsetter.shortlocs(locations);
            sn = cell(numel(slocs), 1);

            for k = 1:numel(slocs)
                loc = slocs{k};
                sn{k} = [dict.(loc{1}), '_', dict.(loc{2})];
            end

            if join
                sn = strjoin(sn, ' ');
            end

        end

        function vals = setFromVector(paramsetter, X)
        % convert from vector to struct representation of the parameter 
            slocs = paramsetter.shortlocs(paramsetter.locations());
            vals = struct();

            for k = 1:numel(slocs)
                s = slocs{k};
                vals = setfield(vals, s{:}, X(k));
            end

        end

        function X = setToVector(paramsetter, vals)
        % convert from struct to vector representation of the parameter 
            locs = paramsetter.shortlocs(paramsetter.locations());
            X = nan(numel(locs), 1);

            for k = 1:numel(locs)
                loc = locs{k};
                X(k) = getfield(vals, loc{:});
            end

        end

        function simsetup = setValues(paramsetter, simsetup, X)
        % Given a parameter vector X, change the model parameter values accordingly.

            model = simsetup.model;
            
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            elyte   = 'Electrolyte';
            co      = 'Coating';
            am      = 'ActiveMaterial';
            sd      = 'SolidDiffusion';
            itf     = 'Interface';

            % Automatic update
            locs = paramsetter.locations();

            for k = 1:numel(locs)
                loc = locs{k};
                model = setfield(model, loc{:}, X(k));
            end

            sn = paramsetter.shortnames;
            
            % Update dependencies (manual work)
            eldes    = {ne, pe};
            sn_eldes = {'ne', 'pe'};
            
            for ielde = 1 : numel(eldes)

                elde    = eldes{ielde};
                sn_elde = sn_eldes{ielde};
                
                if ismember(sprintf('%s_bg', sn_elde), sn)
                    bg    = model.(elde).(co).bruggemanCoefficient;
                    kappa = model.(elde).(co).electronicConductivity;
                    vf    = model.(elde).(co).volumeFraction;
                    model.(elde).(co).effectiveElectronicConductivity = kappa*vf^bg;
                end
                
                if ismember(sprintf('%s_vsa', elde), sn)
                    model.(elde).(co).(am).(sd).volumetricSurfaceArea = model.(elde).(co).(am).(itf).volumetricSurfaceArea;
                end
                

            end

            if model.(elyte).useRegionBruggemanCoefficients

                error('does not support yet active index')
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

            simsetup.model = model;

        end

        function X = getValues(paramsetter, simsetup)
        % Fetch in the model the parameter values and return those as a vector

            model = simsetup.model;
            
            locs  = paramsetter.locations();
            short = paramsetter.shortlocs(locs);

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
