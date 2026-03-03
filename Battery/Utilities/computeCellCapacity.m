function [capacity, capacities] = computeCellCapacity(model)
%
%
% SYNOPSIS:
%   function [capacity, capacities] = computeCellCapacity(model, varargin)
%
% DESCRIPTION: computes the cell usable capacity in Coulomb
%
% PARAMETERS:
%   model - battery model
%
% RETURNS:
%   capacity   - capacity
%   capacities - struct with fields
%      - negativeElectrode
%      - positiveElectrode
%
% EXAMPLE:
%
% SEE ALSO:
%
    ne  = 'NegativeElectrode';
    pe  = 'PositiveElectrode';
    co  = 'Coating';
    itf = 'Interface';
    sd  = 'SolidDiffusion';

    eldes = {ne, pe};

    for ielde = 1 : numel(eldes)

        elde = eldes{ielde};

        if model.(elde).(co).activeMaterialModelSetup.composite

            am1 = 'ActiveMaterial1';
            am2 = 'ActiveMaterial2';

            ams = {am1, am2};

            cap_usable{ielde} = 0;

            vol_fraction = model.(elde).(co).volumeFraction;

            for iam = 1 : numel(ams)

                amc = ams{iam};

                itfmodel = model.(elde).(co).(amc).(itf);

                n    = itfmodel.numberOfElectronsTransferred;
                F    = itfmodel.constants.F;
                G    = itfmodel.G;
                cMax = itfmodel.saturationConcentration;

                switch elde
                  case 'NegativeElectrode'
                    thetaMax = itfmodel.guestStoichiometry100;
                    thetaMin = itfmodel.guestStoichiometry0;
                  case 'PositiveElectrode'
                    thetaMax = itfmodel.guestStoichiometry0;
                    thetaMin = itfmodel.guestStoichiometry100;
                  otherwise
                    error('Electrode not recognized');
                end

                amind = model.(elde).(co).compInds.(amc);
                am_fraction  = model.(elde).(co).volumeFractions(amind);

                vol = sum(am_fraction*vol_fraction.*model.(elde).(co).G.getVolumes());

                cap_usable{ielde} = cap_usable{ielde} + (thetaMax - thetaMin)*cMax*vol*n*F;

            end

        else

            am  = 'ActiveMaterial';

            itfmodel = model.(elde).(co).(am).(itf);

            n    = itfmodel.numberOfElectronsTransferred;
            F    = itfmodel.constants.F;
            G    = itfmodel.G;
            cMax = itfmodel.saturationConcentration;

            switch elde
                
              case 'NegativeElectrode'
                thetaMax = itfmodel.guestStoichiometry100;
                thetaMin = itfmodel.guestStoichiometry0;
              case 'PositiveElectrode'
                thetaMax = itfmodel.guestStoichiometry0;
                thetaMin = itfmodel.guestStoichiometry100;
              otherwise
                error('Electrode not recognized');
            end

            if getJsonStructField(model.(elde).(co).activeMaterialModelSetup, 'swelling', false)
                % Special setup in case of swelling material
                % the guest stochiometries are interpretated as fill-in levels
                
                cMaxTot = model.(elde).(co).maximumTotalConcentration;
                vol = sum(model.(elde).(co).G.getVolumes());
                
                cap_usable{ielde} = (thetaMax - thetaMin)*cMaxTot*vol*n*F;
                
            else
                
                vol_fraction = model.(elde).(co).volumeFraction;

                amind = model.(elde).(co).compInds.(am);
                am_fraction  = model.(elde).(co).volumeFractions(amind);

                vol = sum(am_fraction*vol_fraction.*model.(elde).(co).G.getVolumes());

                cap_usable{ielde} = (thetaMax - thetaMin)*cMax*vol*n*F;

            end
        end


    end

    capacities.(ne) = cap_usable{1};
    capacities.(pe) = cap_usable{2};

    capacity = capacities.(ne);
    ind = value(capacities.(ne)) >= value(capacities.(pe));
    if any(ind)
        capacity(ind) = capacities.(pe)(ind);
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
