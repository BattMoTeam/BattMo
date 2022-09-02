function initstate = initStateChen2020(model, c_ne, c_pe)

    % Some shorthands used for the sub-models
    ne    = 'NegativeElectrode';
    pe    = 'PositiveElectrode';
    am    = 'ActiveMaterial';
    sd    = 'SolidDiffusion';
    elyte = 'Electrolyte';
    itf   = 'Interface';
    ctrl  = 'Control';

    nc = model.G.cells.num;
    T = model.initT;
    initstate.ThermalModel.T = T*ones(nc, 1);

    bat = model;

    initstate = model.updateTemperature(initstate);

    % we setup negative electrode initial state
    
    nitf = bat.(ne).(am).(itf); 

    % We bypass the solid diffusion equation to set directly the particle surface concentration
    switch model.(ne).(am).diffusionModelType
      case 'simple'
        nenp = model.(ne).(am).G.cells.num;
        initstate.(ne).(am).c = c_ne*ones(nenp, 1);
      case 'full'
        nenr = model.(ne).(am).(sd).N;
        nenp = model.(ne).(am).(sd).np;
        initstate.(ne).(am).(sd).c = c_ne*ones(nenr*nenp, 1);
      otherwise
        error('diffusionModelType type not recognized');
    end
    initstate.(ne).(am).(sd).cSurface = c_ne*ones(nenp, 1);

    initstate.(ne).(am) = model.(ne).(am).updateConcentrations(initstate.(ne).(am));
    initstate.(ne).(am).(itf) = nitf.updateOCP(initstate.(ne).(am).(itf));

    OCP = initstate.(ne).(am).(itf).OCP;
    ref = OCP(1);

    initstate.(ne).(am).phi = OCP - ref;

    % we setup positive electrode initial state

    pitf = bat.(pe).(am).(itf); 

    switch model.(ne).(am).diffusionModelType
      case 'simple'
        penp = model.(pe).(am).G.cells.num;
        initstate.(pe).(am).c = c_pe*ones(penp, 1);
      case 'full'
        penr = model.(pe).(am).(sd).N;
        penp = model.(pe).(am).(sd).np;
        initstate.(pe).(am).(sd).c = c_pe*ones(penr*penp, 1);
      otherwise
        error('diffusionModelType type not recognized');
    end

    initstate.(pe).(am).(sd).cSurface = c_pe*ones(penp, 1);
        
    initstate.(pe).(am) = model.(pe).(am).updateConcentrations(initstate.(pe).(am));
    initstate.(pe).(am).(itf) = pitf.updateOCP(initstate.(pe).(am).(itf));

    OCP = initstate.(pe).(am).(itf).OCP;

    initstate.(pe).(am).phi = OCP - ref;

    initstate.(elyte).phi = zeros(bat.(elyte).G.cells.num, 1) - ref;
    initstate.(elyte).c = 1000*ones(bat.(elyte).G.cells.num, 1);

    % setup initial positive electrode external coupling values

    initstate.(ctrl).E = OCP(1) - ref;
    initstate.(ctrl).I = 0;
    initstate.(ctrl).ctrlType = 'constantCurrent';

end
