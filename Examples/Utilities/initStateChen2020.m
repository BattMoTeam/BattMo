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

    nc = model.G.getNumberOfCells();
    T = model.initT;
    initstate.(thermal).T = T*ones(nc, 1);

    bat = model;
    initstate = model.updateTemperature(initstate);

    % we setup negative electrode initial state

    nitf = bat.(ne).(co).(am).(itf);

    % We bypass the solid diffusion equation to set directly the particle surface concentration
    switch model.(ne).(co).(am).diffusionModelType
      case 'simple'
        nenp = model.(ne).(co).G.getNumberOfCells();
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
        penp = model.(pe).(co).G.getNumberOfCells();
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

    initstate.(elyte).phi = zeros(bat.(elyte).G.getNumberOfCells(), 1) - ref;
    initstate.(elyte).c = 1000*ones(bat.(elyte).G.getNumberOfCells(), 1);

    % setup initial positive electrode external coupling values

    initstate.(ctrl).E = OCP(1) - ref;
    initstate.(ctrl).I = 0;
    initstate.(ctrl).ctrlType = 'constantCurrent';

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
