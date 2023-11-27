function initstate = initStateChen2020(model, c_ne, c_pe)

    % Some shorthands used for the sub-models
    ne      = 'NegativeElectrode';
    pe      = 'PositiveElectrode';
    am      = 'ActiveMaterial';
    sd      = 'SolidDiffusion';
    elyte   = 'Electrolyte';
    itf     = 'Interface';
    ctrl    = 'Control';
    co      = 'Coating';
    thermal = 'ThermalModel';

    nc = model.G.cells.num;
    T = model.initT;
    initstate.(thermal).T = T*ones(nc, 1);

    bat = model;
    initstate = model.updateTemperature(initstate);

    % we setup negative electrode initial state

    nitf = bat.(ne).(co).(am).(itf);

    % We bypass the solid diffusion equation to set directly the particle surface concentration
    switch model.(ne).(co).(am).diffusionModelType
      case 'simple'
        nenp = model.(ne).(co).G.cells.num;
        initstate.(ne).(co).(am).(sd).cAverage = c_ne*ones(nenp, 1);
        initstate.(ne).(co).(am).(sd).cSurface = c_ne*ones(nenp, 1);
      case 'full'
        nenr = model.(ne).(co).(am).(sd).N;
        nenp = model.(ne).(co).(am).(sd).np;
        initstate.(ne).(co).(am).(sd).c = c_ne*ones(nenr*nenp, 1);
      otherwise
        error('diffusionModelType type not recognized');
    end
    initstate.(ne).(co).(am).(sd).cSurface = c_ne*ones(nenp, 1);

    initstate.(ne).(co).(am) = model.(ne).(co).(am).updateConcentrations(initstate.(ne).(co).(am));

    initstate = model.evalVarName(initstate, {ne, co, am, itf, 'OCP'});
    initstate.(ne).(co).(am).(itf) = nitf.updateOCP(initstate.(ne).(co).(am).(itf));
    OCP = initstate.(ne).(co).(am).(itf).OCP;
    ref = OCP(1);

    initstate.(ne).(co).phi = OCP - ref;

    % we setup positive electrode initial state

    pitf = bat.(pe).(co).(am).(itf);

    switch model.(ne).(co).(am).diffusionModelType
      case 'simple'
        penp = model.(pe).(co).G.cells.num;
        initstate.(pe).(co).(am).(sd).cAverage = c_pe*ones(penp, 1);
        initstate.(pe).(co).(am).(sd).cSurface = c_pe*ones(penp, 1);
      case 'full'
        penr = model.(pe).(co).(am).(sd).N;
        penp = model.(pe).(co).(am).(sd).np;
        initstate.(pe).(co).(am).(sd).c = c_pe*ones(penr*penp, 1);
      otherwise
        error('diffusionModelType type not recognized');
    end

    initstate.(pe).(co).(am).(sd).cSurface = c_pe*ones(penp, 1);

    initstate.(pe).(co).(am) = model.(pe).(co).(am).updateConcentrations(initstate.(pe).(co).(am));

    initstate = model.evalVarName(initstate, {pe, co, am, itf, 'OCP'});
    initstate.(pe).(co).(am).(itf) = pitf.updateOCP(initstate.(pe).(co).(am).(itf));
    OCP = initstate.(pe).(co).(am).(itf).OCP;
    initstate.(pe).(co).phi = OCP - ref;

    initstate.(elyte).phi = zeros(bat.(elyte).G.cells.num, 1) - ref;
    initstate.(elyte).c = 1000*ones(bat.(elyte).G.cells.num, 1);

    % setup initial positive electrode external coupling values

    initstate.(ctrl).E = OCP(1) - ref;
    initstate.(ctrl).I = 0;
    initstate.(ctrl).ctrlType = 'constantCurrent';

end
