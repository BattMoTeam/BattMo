function [model, initstate] = setupPlatingInitialState(model, T, cElectrolyte, phiElectrolyte, cElectrodeInit, Imax)

    %%
    % We define some shortcuts

    ne      = 'NegativeElectrode';
    pe      = 'PositiveElectrode';
    lp      = 'LithiumPlating';
    elyte   = 'Electrolyte';
    thermal = 'ThermalModel';
    co      = 'Coating';
    am      = 'ActiveMaterial';
    itf     = 'Interface';
    sd      = 'SolidDiffusion';
    ctrl    = 'Control';
    cc      = 'CurrentCollector';
    
    N = model.(sd).N;
    initstate.(sd).c        = cElectrodeInit*ones(N, 1);
    initstate.(sd).cSurface = cElectrodeInit;

    initstate.T = T;
    initstate.(itf).cElectrolyte   = cElectrolyte;
    initstate.(itf).phiElectrolyte = phiElectrolyte;

    initstate = model.evalVarName(initstate, {itf, 'OCP'});
    OCP = initstate.(itf).OCP;
    initstate.E = OCP + phiElectrolyte;

    F = model.(itf).constants.F;
    R = model.(itf).constants.R;

    if model.useLithiumPlating
        
        thresholdParameter   = model.LithiumPlating.thresholdParameter;
        r                    = model.LithiumPlating.particleRadius;
        vf                   = model.LithiumPlating.volumeFraction;
        platedConcentration0 = thresholdParameter * vf / ((4/3)*pi*r^3);
        
        %%
        % initialisation so that the overpotential are zero at the beginning
        platedConcentrationInit = platedConcentration0/(exp((F*OCP)/(R*T)) - 1)^(1/4);

        model.(lp).platedReferenceConcentration = platedConcentrationInit;

        initstate.(lp).platedConcentration     = platedConcentrationInit ;
        initstate.(lp).platedConcentrationNorm = platedConcentrationInit/model.(lp).platedReferenceConcentration;
        initstate.(lp).phiSolid                = initstate.E;
        initstate.(lp).phiElectrolyte          = phiElectrolyte;
        initstate.(lp).cElectrolyte            = cElectrolyte;
        initstate.(lp).nSEI                    = 0;
        
    end

    scalingparams = struct('I'                  , Imax                              , ...
                           'elyteConcentration' , initstate.(itf).cElectrolyte);

    if model.useLithiumPlating
        scalingparams.platedConcentration = model.(lp).platedReferenceConcentration;
    end

    %%
    % We setup the scaling for the residual equations
    %

    model = model.setupScalings(scalingparams);
    
end

%{
Copyright 2021-2026 SINTEF Industry, Sustainable Energy Technology
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
