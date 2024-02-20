classdef NPlengthSetter1D

    properties

        NPratio % NP ratio

        lengthSetter
        porositySetter

        % helpers
        cRatio   % For a given NP ratio, cRatio  = (volumefraction(ne)*length(ne))/(volumefraction(pe)*length(pe)) is a constant
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
            nplengthsetter.cRatio         = nplengthsetter.computeMultiplicationCoefficient(model);
            nplengthsetter.lengthSetter   = lengthSetter;
            nplengthsetter.porositySetter = porositySetter;

        end

        function model = setValues(nplengthsetter, model, values)

            ctrl = 'Control';

            fdinds         = nplengthsetter.fdinds;
            cRatio         = nplengthsetter.cRatio;
            lengthSetter   = nplengthsetter.lengthSetter;
            porositySetter = nplengthsetter.porositySetter;

            pelength = values(fdinds.pelength);
            neporo   = values(fdinds.neporo);
            peporo   = values(fdinds.peporo);

            v = double2ADI(zeros(2, 1), pelength);
            v(1) = cRatio*(1 - peporo)./(1 - neporo).*pelength;
            v(2) = pelength;

            model = lengthSetter.setLengths(model, v);
            model = porositySetter.setValues(model, [neporo; peporo]);

            C = computeCellCapacity(model);

            DRate = model.(ctrl).DRate;
            model.(ctrl).Imax = (C/hour)*DRate;

        end

        function values = getValues(nplengthsetter, model)

            lengthSetter   = nplengthsetter.lengthSetter;
            porositySetter = nplengthsetter.porositySetter;

            pelength = lengthSetter.getAllLengths(model);
            pelength = pelength(3);

            poros = porositySetter.getValues(model);

            values = [pelength; poros];

        end


        function cRatio = computeMultiplicationCoefficient(nplengthsetter, model)

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
                idx      = model.(elde).(co).compInds.(am);
                amvf     = model.(elde).(co).volumeFractions(idx);

                F        = model.(elde).(co).constants.F;

                d.(elde) = abs(theta100 - theta0)*cmax*amvf*F/hour;

            end

            cRatio = d.(pe)/d.(ne)*NPratio;

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
