function [cap, cap_neg, cap_pos, specificEnergy] = computeCellCapacity(model, varargin)

%
%
% SYNOPSIS:
%   function c = computeCellCapacity(model, varargin)
%
% DESCRIPTION: computes the cell usable capacity in Coulomb
%
% PARAMETERS:
%   model - battery model
%
% RETURNS:
%   cap - capacity
%
% EXAMPLE:
%
% SEE ALSO:
%
    opt = struct('packingMass', 0);
    opt = merge_options(opt, varargin{:});

    ne  = 'NegativeElectrode';
    pe  = 'PositiveElectrode';
    am  = 'ActiveMaterial';
    itf = 'Interface';
    sd  = 'SolidDiffusion';
    
    eldes = {ne, pe};
    
    for ind = 1 : numel(eldes)
        
        elde = eldes{ind};
        
        ammodel = model.(elde).(am);
        itfmodel = ammodel.(itf);
        
        n    = itfmodel.n;
        F    = itfmodel.constants.F;
        G    = itfmodel.G;
        cMax = itfmodel.cmax;

        switch elde
          case 'NegativeElectrode'
            thetaMax = itfmodel.theta100;
            thetaMin = itfmodel.theta0;
          case 'PositiveElectrode'
            thetaMax = itfmodel.theta0;
            thetaMin = itfmodel.theta100;            
          otherwise
            error('Electrode not recognized');
        end
        
        vol_fraction = ammodel.volumeFraction;
        am_fraction  = ammodel.activeMaterialFraction;
        
        vol = sum(am_fraction*vol_fraction.*ammodel.G.cells.volumes);
        
        cap_usable(ind) = (thetaMax - thetaMin)*cMax*vol*n*F;
        
    end
    
    cap_neg = cap_usable(1);
    cap_pos = cap_usable(2);
    
    cap = min(cap_usable); 

    
    if nargout > 3
        
        r = cap_neg/cap_pos;
        
        thetaMinPos = model.(pe).(am).(itf).theta100;
        thetaMaxPos = model.(pe).(am).(itf).theta0;
        thetaMinNeg = model.(ne).(am).(itf).theta0;
        thetaMaxNeg = model.(ne).(am).(itf).theta100;
        
        elde = 'PositiveElectrode';

        ammodel = model.(elde).(am);
        itfmodel = model.(elde).(am).(itf);
        
        F = itfmodel.constants.F;
        G = itfmodel.G;
        n = itfmodel.n;
        assert(n == 1, 'not implemented yet');
        cMax = itfmodel.cmax;
        
        vol_fraction = ammodel.volumeFraction;
        am_fraction  = ammodel.activeMaterialFraction;
        vol = sum(am_fraction*vol_fraction.*ammodel.G.cells.volumes);
        
        func = @(theta) model.(elde).(am).(itf).computeOCPFunc(theta, 298, 1);

        thetaMax = min(thetaMaxPos, thetaMinPos + r*(thetaMaxPos - thetaMinPos));

        theta = linspace(thetaMinPos, thetaMax, 100);
        energy = sum(func(theta(1 : end - 1)).*diff(theta)*vol*F*cMax);
        
        elde = 'NegativeElectrode';        
        
        itfmodel = model.(elde).(am).(itf);
        F = itfmodel.constants.F;
        G = itfmodel.G;
        n = itfmodel.n;
        assert(n == 1, 'not implemented yet');
        cMax = itfmodel.cmax;
        
        vol_fraction = ammodel.volumeFraction;
        am_fraction  = ammodel.activeMaterialFraction;
        vol = sum(am_fraction*vol_fraction.*ammodel.G.cells.volumes);
        
        func = @(theta) model.(elde).(am).(itf).computeOCPFunc(theta, 298, 1);

        thetaMin = max(thetaMinNeg, thetaMaxNeg - 1/r*(thetaMaxNeg - thetaMinNeg));

        theta = linspace(thetaMin, thetaMaxNeg, 100);
        energy = energy - sum(func(theta(1 : end - 1)).*diff(theta)*vol*F*cMax);
        
        mass = computeCellMass(model, 'packingMass', opt.packingMass);
        
        specificEnergy = energy/mass;
        
    else
        
        specificEnergy = [];
        
    end
    
end



%{
Copyright 2021-2022 SINTEF Industry, Sustainable Energy Technology
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
