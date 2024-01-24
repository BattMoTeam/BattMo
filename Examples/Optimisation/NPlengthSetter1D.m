classdef NPlengthSetter1D

    properties

        NPratio % NP ratio

        lengthSetter
        porositySetter

        % helpers
        alpha   % For a given NP ratio, alpha is equal the corresponding ratio between the thicknesses of the electrodes
        fdnames
        fdinds

    end

    methods

        function nplengthsetter = NPlengthSetter1D(model, gridGenerator, NPratio)

            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';

            lengthSetter   = LengthSetter1D(gridGenerator, {ne, pe});
            porositySetter = PorositySetter(model, {ne, pe});

            fdnames = {'pelength', 'neporo', 'peporo'};
            for ifd = 1 : numel(fdnames)
                fdname = fdnames{ifd};
                fdinds.(fdname) = ifd;
            end

            nplengthsetter.NPratio        = NPratio;
            nplengthsetter.fdinds         = fdinds;
            nplengthsetter.fdnames        = fdnames;
            nplengthsetter.alpha          = nplengthsetter.computeMultiplicationCoefficient(model);
            nplengthsetter.lengthSetter   = lengthSetter;
            nplengthsetter.porositySetter = porositySetter;

        end

        function model = setValues(nplengthsetter, model, values)

            ctrl = 'Control';

            fdinds         = nplengthsetter.fdinds;
            alpha          = nplengthsetter.alpha;
            lengthSetter   = nplengthsetter.lengthSetter;
            porositySetter = nplengthsetter.porositySetter;

            pelength = values(fdinds.pelength);
            neporo   = values(fdinds.neporo);
            peporo   = values(fdinds.peporo);

            v = double2ADI(zeros(2, 1), pelength);
            v(1) = alpha*(1 - peporo)./(1 - neporo).*pelength;
            v(2) = pelength;

            model = lengthSetter.setLengths(model, v);
            model = porositySetter.setValues(model, [neporo; peporo]);

            C = computeCellCapacity(model);

            CRate = model.(ctrl).CRate;
            model.(ctrl).Imax = (C/hour)*CRate;

        end

        function values = getValues(nplengthsetter, model)

            lengthSetter   = nplengthsetter.lengthSetter;
            porositySetter = nplengthsetter.porositySetter;

            pelength = lengthSetter.getAllLengths(model);
            pelength = pelength(3);

            poros = porositySetter.getValues(model);

            values = [pelength; poros];

        end


        function alpha = computeMultiplicationCoefficient(nplengthsetter, model)

            NPratio = nplengthsetter.NPratio;

            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            am  = 'ActiveMaterial';
            itf = 'Interface';
            co  = 'Coating';

            eldes = {ne, pe};

            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};

                theta100 = model.(elde).(co).(am).(itf).guestStoichiometry100;
                theta0   = model.(elde).(co).(am).(itf).guestStoichiometry0;
                cmax     = model.(elde).(co).(am).(itf).saturationConcentration;
                eldevf   = model.(elde).(co).volumeFraction;
                idx      = model.(elde).(co).compInds.(am);
                amvf     = model.(elde).(co).volumeFractions(idx);

                F        = model.(elde).(co).constants.F;

                d.(elde) = abs(theta100 - theta0)*cmax*eldevf*amvf*F/hour;

            end

            alpha = d.(pe)/d.(ne)*NPratio;

        end
    end

end
