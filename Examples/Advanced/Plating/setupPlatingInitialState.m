function [model, initstate] = setupPlatingInitialState(model, T, cElectrolyte, phiElectrolyte, cElectrodeInit)

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
        
        nPl0                 = model.LithiumPlating.nPl0;
        r                    = model.LithiumPlating.particleRadius;
        vf                   = model.LithiumPlating.volumeFraction;
        platedConcentration0 = nPl0 * vf / ((4/3)*pi*r^3);
        
        %%
        % initialisation so that the overpotential are zero at the beginning
        platedConcentrationInit = platedConcentration0/(exp((F*OCP)/(R*T)) - 1)^(1/4);

        model.(lp).platedConcentrationRef = platedConcentrationInit;

        initstate.(lp).platedConcentration     = platedConcentrationInit ;
        initstate.(lp).platedConcentrationNorm = platedConcentrationInit/model.(lp).platedConcentrationRef;
        initstate.(lp).phiSolid                = initstate.E;
        initstate.(lp).phiElectrolyte          = phiElectrolyte;
        initstate.(lp).cElectrolyte            = cElectrolyte;
        initstate.(lp).nSEI                    = 0;
        
    end

end
