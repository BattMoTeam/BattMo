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
    am  = 'ActiveMaterial';
    itf = 'Interface';
    sd  = 'SolidDiffusion';
    
    eldes = {ne, pe};
    
    for ind = 1 : numel(eldes)
        
        elde = eldes{ind};

        switch model.(elde).(co).activematerial_type

          case 'default'
            
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
            
            vol_fraction = model.(elde).(co).volumeFraction;
            
            amind = model.(elde).(co).compInds.(am);
            am_fraction  = model.(elde).(co).volumeFractions(amind);
            
            vol = sum(am_fraction*vol_fraction.*model.(elde).(co).G.cells.volumes);
            
            cap_usable{ind} = (thetaMax - thetaMin)*cMax*vol*n*F;
            
          case 'composite'

            % we know we deal with a negative electrode
            ammodel = model.(elde).(am);

            gr = 'FirstMaterial';
            si = 'SecondMaterial';

            mats = {gr, si};

            cap_usable{ind} = 0;
            
            for imat = 1 : numel(mats)
                
                mat = mats{imat};

                itfmodel = ammodel.(mat).(itf);
                n    = itfmodel.n;
                F    = itfmodel.constants.F;
                G    = itfmodel.G;
                cMax = itfmodel.cmax;

                thetaMax = itfmodel.guestStoichiometry100;
                thetaMin = itfmodel.guestStoichiometry0;

                vol_fraction = ammodel.volumeFraction;
                am_fraction  = ammodel.(mat).activeMaterialFraction;
                
                vol = sum(am_fraction*vol_fraction.*ammodel.G.cells.volumes);

                cap_usable{ind} = cap_usable{ind} + (thetaMax - thetaMin)*cMax*vol*n*F;
            end
            
          otherwise
            error('electrode_case not recognized');
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
