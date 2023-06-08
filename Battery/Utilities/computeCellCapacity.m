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
    opt = struct('packingMass', 0, ...
                'isSwellingMaterial', false);
    opt = merge_options(opt, varargin{:});

    ne  = 'NegativeElectrode';
    pe  = 'PositiveElectrode';
    am  = 'ActiveMaterial';
    itf = 'Interface';
    sd  = 'SolidDiffusion';
    
    eldes = {ne, pe};
    
    for ind = 1 : numel(eldes)
        
        elde = eldes{ind};

        switch model.(elde).electrode_case

          case 'default'
            
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
            
          case 'composite'

            % we know we deal with a negative electrode
            ammodel = model.(elde).(am);

            gr = 'FirstMaterial';
            si = 'SecondMaterial';

            mats = {gr, si};

            cap_usable(ind) = 0;
            
            for imat = 1 : numel(mats)
                
                mat = mats{imat};

                itfmodel = ammodel.(mat).(itf);
                n    = itfmodel.n;
                F    = itfmodel.constants.F;
                G    = itfmodel.G;
                cMax = itfmodel.cmax;

                thetaMax = itfmodel.theta100;
                thetaMin = itfmodel.theta0;

                vol_fraction = ammodel.volumeFraction;
                am_fraction  = ammodel.(mat).activeMaterialFraction;
                
                vol = sum(am_fraction*vol_fraction.*ammodel.G.cells.volumes);

                cap_usable(ind) = cap_usable(ind) + (thetaMax - thetaMin)*cMax*vol*n*F;
            end
            
          otherwise
            error('electrode_case not recognized');
        end
        
        % Warning : The calculation we did above relies on the initialstate
        % of the activeMaterial (the one defined in the json file).
        % However, for a swelling material, the calculation of the capacity
        % is a bit more complicated : if the original soc is not 0, the
        % particle is already a bit lithiated and thus has already swelled. Therefore, it has no sense to say
        % that they can host cmx Li atoms; however, it is still true that
        % they can host Nmax = cmax * 4/3*pi*R_delith^3. Still, we have
        % that the real capacity is Nmax_part * NtotParticles with NtotParticles = totalVolume/particleVolume. 


        if opt.isSwellingMaterial && ammodel.isSwellingMaterial
            vol = sum(vol_fraction.*ammodel.G.cells.volumes);

            R_delith = ammodel.SolidDiffusion.rp;
            soc = 0.9;
           
            %to get the radius of the particle unlithiated, we tuse the eq
            %11 in Chandrasekaran, Magasinski, Yushin and Fuller
            molarVolumeSi = 1.2e-05;
            molarVolumeLi = 9e-06;
            Q = (3.75.*molarVolumeLi)./(molarVolumeSi);

            R_initial = R_delith * (1 + Q * soc)^(1/3);

            Nmax_part = cMax * (4/3) * pi * (1+Q) * R_delith^3;

            Npart_tot = vol / ((4/3) * pi * R_initial^3);

            cap_usable(ind) =  Npart_tot *  Nmax_part * n * F ;
        end
        
    end
 
    cap_neg = cap_usable(1)
    cap_pos = cap_usable(2)

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

        theta = linspace(thetaMinPos, thetaMax, 1000);
        energy = sum(func(theta(1 : end - 1)).*diff(theta)*vol*F*cMax);
        
        elde = 'NegativeElectrode';        

        ammodel  = model.(elde).(am);
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

        theta = linspace(thetaMin, thetaMaxNeg, 1000);

        energy = energy - sum(func(theta(1 : end - 1)).*diff(theta)*vol*F*cMax);
        
        mass = computeCellMass(model, 'packingMass', opt.packingMass);
        
        specificEnergy = energy/mass;
        
    else
        
        specificEnergy = [];
        
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
