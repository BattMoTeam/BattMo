classdef PorositySetter

    properties

        % helper properties

        compnames % names of the components
        compinds  % index of the components
        elytemaps
        refvalues

    end

    methods

        function porosetter = PorositySetter(model, components)

            ne    = 'NegativeElectrode';
            sep   = 'Separator';
            pe    = 'PositiveElectrode';
            elyte = 'Electrolyte';
            co    = 'Coating';

            if model.use_thermal
                error('Porisity setter is not implemented for thermal model yet (not difficult but not done...)')
            end

            compnames = {ne, sep, pe};

            compinds = ismember({ne, sep, pe}, components);

            elyteinds = zeros(model.(elyte).G.parentGrid.getNumberOfCells(), 1);
            elyteinds(model.(elyte).G.mappings.cellmap) = (1 : model.(elyte).G.getNumberOfCells())';

            refvalues = zeros(numel(compnames), 1);

            for icomp = 1 : numel(compnames)
                compname = compnames{icomp};
                switch compname
                  case ne
                    compmodel = model.(ne).(co);
                    poro = 1 - compmodel.volumeFraction;
                  case sep
                    compmodel = model.(sep);
                    poro = compmodel.porosity;
                  case pe
                    compmodel = model.(pe).(co);
                    poro = 1 - compmodel.volumeFraction;
                  otherwise
                    error('compname not recognized');
                end
                refvalues(icomp) = poro;
                elytemaps.(compname) = elyteinds(compmodel.G.mappings.cellmap);

            end

            porosetter.compnames = compnames;
            porosetter.compinds  = compinds;
            porosetter.refvalues = refvalues;
            porosetter.elytemaps = elytemaps;

        end

        function model = setPorosities(porosetter, model, v)

            compnames = porosetter.compnames;
            compinds  = porosetter.compinds;
            refvalues = porosetter.refvalues;

            poros = refvalues;
            poros = subsasgnAD(poros, compinds, v);

            for icomp = 1 : numel(compnames)
                compname = compnames{icomp};
                model = porosetter.setComponentPorosity(model, poros(icomp), compname);
            end

        end

        function model = setComponentPorosity(porosetter, model, poro, compname)
        % We update the porosities and the dependent values.

            elytemaps = porosetter.elytemaps;

            ne    = 'NegativeElectrode';
            sep   = 'Separator';
            pe    = 'PositiveElectrode';
            elyte = 'Electrolyte';
            co    = 'Coating';
            am    = 'ActiveMaterial';
            itf   = 'Interface';
            sd    = 'SolidDiffusion';

            switch compname

              case {ne, pe}

                model.(compname).(co).volumeFraction = 1 - poro;
                switch model.(compname).(co).(am).diffusionModelType
                  case {'simple', 'interParticleOnly'}
                    % nothing to do
                  case 'full'
                    model.(compname).(co).(am).(sd).volumeFraction = 1 - poro;
                  otherwise
                    error('Unknown diffusionModelType %s', diffusionModelType);
                end

                % We need to update the effective conductivities

                kappa = model.(compname).(co).electronicConductivity;
                vf    = model.(compname).(co).volumeFraction;
                bg    = model.(compname).(co).bruggemanCoefficient;

                model.(compname).(co).effectiveElectronicConductivity = kappa*vf^bg;

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

        function v = getPorosities(porosetter, model)

            compinds = porosetter.compinds;

            v = porosetter.getAllPorosities(model);
            v = v(compinds);

        end
    end

end
