classdef NPlengthSetter1D

    properties

        NPratio % NP ratio

        lengthSetter

        % helpers
        cRatio   % For a given NP ratio, cRatio  = (volumefraction(ne)*length(ne))/(volumefraction(pe)*length(pe)) is a constant
        volumeFractions

    end

    methods

        function nplengthsetter = NPlengthSetter1D(model, gridGenerator, NPratio)

            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            co  = 'Coating';
            
            lengthSetter   = LengthSetter1D(gridGenerator, {ne, pe});


            volumeFractions.(ne) = model.(ne).(co).volumeFraction;
            volumeFractions.(pe) = model.(pe).(co).volumeFraction;
            
            nplengthsetter.NPratio         = NPratio;
            nplengthsetter.cRatio          = nplengthsetter.computeMultiplicationCoefficient(model);
            nplengthsetter.lengthSetter    = lengthSetter;
            nplengthsetter.volumeFractions = volumeFractions;
            
        end

        function model = setValue(nplengthsetter, model, value)

            ctrl = 'Control';

            cRatio       = nplengthsetter.cRatio;
            vfs          = nplengthsetter.volumeFractions;
            lengthSetter = nplengthsetter.lengthSetter;
            
            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';

            pelength = value;

            v = double2ADI(zeros(2, 1), pelength);

            v(1) = cRatio*(vfs.(pe)./vfs.(ne)).*pelength;
            v(2) = pelength;

            model = lengthSetter.setLengths(model, v);

            C = computeCellCapacity(model);

            DRate = model.(ctrl).DRate;
            model.(ctrl).Imax = (C/hour)*DRate;

        end

        function value = getValue(nplengthsetter, model)

            lengthSetter   = nplengthsetter.lengthSetter;

            pelength = lengthSetter.getAllLengths(model);
            pelength = pelength(3);

            value = pelength;

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
