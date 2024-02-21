classdef PorositySetter

    properties

        % helper properties

        compnames % names of the components
        compinds  % index of the components
        elytemaps
        refvalues
        
        sumspecvols % Sum of the specific volumes for each electrodes. Used to set effective density of the coating.
        
    end

    methods

        function porosetter = PorositySetter(model, components)

            ne    = 'NegativeElectrode';
            sep   = 'Separator';
            pe    = 'PositiveElectrode';
            elyte = 'Electrolyte';
            co    = 'Coating';
            am    = 'ActiveMaterial';
            bd    = 'Binder';
            ad    = 'ConductingAdditive';

            if model.use_thermal
                error('Porisity setter is not implemented for thermal model yet (not difficult but not done...)')
            end

            compnames = {ne, sep, pe};

            compinds = ismember({ne, sep, pe}, components);

            elyteinds = zeros(model.(elyte).G.parentGrid.topology.cells.num, 1);
            elyteinds(model.(elyte).G.mappings.cellmap) = (1 : model.(elyte).G.getNumberOfCells())';

            refvalues = nan(numel(compnames), 1);

            mats = {am, bd, ad};
            
            for icomp = 1 : numel(compnames)
                compname = compnames{icomp};
                switch compname
                  case {ne, pe}
                    compmodel = model.(compname).(co);
                    poro = 1 - compmodel.volumeFraction;
                    sumspecvols.(compname) = 0;
                    for imat = 1 : numel(mats)
                        mat = mats{imat};
                        sumspecvols.(compname) = sumspecvols.(compname) + model.(compname).(co).(mat).massFraction/model.(compname).(co).(mat).density;
                    end
                  case sep
                    compmodel = model.(sep);
                    poro = compmodel.porosity;
                  otherwise
                    error('compname not recognized');
                end
                refvalues(icomp) = poro;
                elytemaps.(compname) = elyteinds(compmodel.G.mappings.cellmap);

            end

            porosetter.compnames   = compnames;
            porosetter.compinds    = compinds;
            porosetter.refvalues   = refvalues;
            porosetter.elytemaps   = elytemaps;
            porosetter.sumspecvols = sumspecvols;
            
        end

        function model = setValues(porosetter, model, v)

            compnames = porosetter.compnames;
            compinds  = porosetter.compinds;
            refvalues = porosetter.refvalues;

            poros = refvalues;
            poros = subsasgnAD(poros, compinds, v);

            for icomp = 1 : numel(compnames)
                compname = compnames{icomp};
                model = porosetter.setupComponent(model, poros(icomp), compname);
            end

        end

        function model = setupComponent(porosetter, model, poro, compname)
        % We update the porosities and the dependent values.
            
            elytemaps = porosetter.elytemaps;

            ne    = 'NegativeElectrode';
            sep   = 'Separator';
            pe    = 'PositiveElectrode';
            elyte = 'Electrolyte';
            co    = 'Coating';
            sd    = 'SolidDiffusion';
            am    = 'ActiveMaterial';
            bd    = 'ActiveMaterial';
            ca    = 'ConductingAdditive';

            switch compname

              case {ne, pe}

                model.(compname).(co).volumeFraction = 1 - poro;

                switch model.(compname).(co).(am).diffusionModelType
                  case 'simple'
                    % nothing to do
                  case 'full'

                    compInds = model.(compname).(co).compInds;
                    amvf     = model.(compname).(co).volumeFractions(compInds.(am));
                    
                    model.(compname).(co).(am).(sd).volumeFraction = (1 - poro)*amvf;
                    
                  otherwise
                    error('Unknown diffusionModelType %s', diffusionModelType);
                end

                % We need to update the effective conductivities

                kappa = model.(compname).(co).electronicConductivity;
                vf    = model.(compname).(co).volumeFraction;
                bg    = model.(compname).(co).bruggemanCoefficient;

                model.(compname).(co).effectiveElectronicConductivity = kappa*vf^bg;

                % we set the effective density of the coating
                model.(compname).(co).effectiveDensity = vf/porosetter.sumspecvols.(compname);
                
              case sep

                model.(sep).porosity = poro;

              otherwise

                error('compname not recognized');

            end

            poro = poro.*ones(numel(elytemaps.(compname)), 1);
            model.(elyte).volumeFraction = subsasgnAD(model.(elyte).volumeFraction, elytemaps.(compname), poro);
            % Note that we do not need to update the electrolyte effective conductivity as it is done dynamically using  state variable.

        end


        function v = getAllPorosities(porosetter, model)

            compnames = porosetter.compnames;

            ne    = 'NegativeElectrode';
            sep   = 'Separator';
            pe    = 'PositiveElectrode';
            elyte = 'Electrolyte';
            co    = 'Coating';
            am    = 'ActiveMaterial';

            v = nan(numel(compnames), 1);

            for icomp = 1 : numel(compnames)

                compname = compnames{icomp};
                switch compname
                  case {ne, pe}
                    poro = 1 - model.(compname).(co).volumeFraction;
                  case sep
                    poro = model.(sep).porosity;
                  otherwise
                    error('compname not recognized');
                end

                v = subsasgnAD(v, icomp, poro(1));

            end

        end

        function v = getValues(porosetter, model)

            compinds = porosetter.compinds;

            v = porosetter.getAllPorosities(model);
            v = v(compinds);

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
