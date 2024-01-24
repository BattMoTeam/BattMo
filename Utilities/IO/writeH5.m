function writeH5(model, states, name)
% WRITEOUTPUT Summary of this function goes here
%   Detailed explanation goes here

    ne      = 'NegativeElectrode';
    pe      = 'PositiveElectrode';
    elyte   = 'Electrolyte';
    thermal = 'ThermalModel';
    co      = 'Coating';
    am      = 'ActiveMaterial';
    itf     = 'Interface';
    sd      = 'SolidDiffusion';
    ctrl    = 'Control';
    cc      = 'CurrentCollector';

    file_name = fullfile(battmoDir(), name);

    % Check if the file exists, and if so, delete it
    if exist(file_name, 'file') == 2
        % Delete the file
        delete(file_name);
    end

    fid = H5F.create(file_name);

    CRate = model.Control.CRate;

    % Instantiate state variables
    time = zeros(size(states));
    voltage = zeros(size(states));
    current = zeros(size(states));

    ne_electricPotential = zeros(length(states), length(states{1}.(ne).(co).phi));
    ne_concentration = zeros(length(states), length(states{1}.(ne).(co).(am).(sd).c));
    ne_surfaceConcentration = zeros(length(states), length(states{1}.(ne).(co).(am).(sd).cSurface));

    pe_electricPotential = zeros(length(states), length(states{1}.(pe).(co).phi));
    pe_concentration = zeros(length(states), length(states{1}.(pe).(co).(am).(sd).c));
    pe_surfaceConcentration = zeros(length(states), length(states{1}.(pe).(co).(am).(sd).cSurface));

    elyte_electricPotential = zeros(length(states), length(states{1}.(elyte).phi));
    elyte_concentration = zeros(length(states), length(states{1}.(elyte).c));

    % Process the states array
    for ind = 1:length(states)
        % Get global state variables
        time(ind, 1) = states{ind}.time;
        voltage(ind, 1) = states{ind}.Control.E;
        current(ind, 1) = states{ind}.Control.I;

        % Get local state variables - Negative Electrode
        ne_electricPotential(ind, :) = states{ind}.(ne).(co).phi;
        ne_surfaceConcentration(ind, :) = states{ind}.(ne).(co).(am).(sd).cSurface;
        ne_concentration(ind, :) = states{ind}.(ne).(co).(am).(sd).c;

        % Get local state variables - Positive Electrode
        pe_electricPotential(ind, :) = states{ind}.(pe).(co).phi;
        pe_surfaceConcentration(ind, :) = states{ind}.(pe).(co).(am).(sd).cSurface;
        pe_concentration(ind, :) = states{ind}.(pe).(co).(am).(sd).c;

        % Get local state variables - Electrolyte
        elyte_electricPotential(ind, :) = states{ind}.(elyte).phi;
        elyte_concentration(ind, :) = states{ind}.(elyte).c;

    end

    % Instantiate global state variables
    h5create(file_name, '/global/quantities/time', size(time));
    h5create(file_name, '/global/quantities/voltage', size(voltage));
    h5create(file_name, '/global/quantities/current', size(current));

    % Instantiate local state variables
    h5create(file_name, '/local/NegativeElectrode/quantities/surfaceConcentration', size(ne_surfaceConcentration));
    h5create(file_name, '/local/NegativeElectrode/quantities/concentration', size(ne_concentration));
    h5create(file_name, '/local/NegativeElectrode/quantities/electricPotential', size(ne_electricPotential));
    h5create(file_name, '/local/PositiveElectrode/quantities/surfaceConcentration', size(pe_surfaceConcentration));
    h5create(file_name, '/local/PositiveElectrode/quantities/concentration', size(pe_concentration));
    h5create(file_name, '/local/PositiveElectrode/quantities/electricPotential', size(pe_electricPotential));
    h5create(file_name, '/local/Electrolyte/quantities/concentration', size(elyte_concentration));
    h5create(file_name, '/local/Electrolyte/quantities/electricPotential', size(elyte_electricPotential));

    % Instantiate specifications
    h5create(file_name, '/global/specification/cRate', size(CRate));

    % Instantiate global grid properties
    h5create(file_name, '/global/grid/cells/num', size(model.G.getNumberOfCells()));
    h5create(file_name, '/global/grid/cells/facePos', size(model.G.cells.facePos));
    h5create(file_name, '/global/grid/cells/indexMap', size(model.G.cells.indexMap));
    h5create(file_name, '/global/grid/cells/faces', size(model.G.cells.faces));
    h5create(file_name, '/global/grid/cells/volumes', size(model.G.getVolumes()));
    h5create(file_name, '/global/grid/cells/centroids', size(model.G.cells.centroids));
    h5create(file_name, '/global/grid/faces/num', size(model.G.getNumberOfFaces()));
    h5create(file_name, '/global/grid/faces/nodePos', size(model.G.faces.nodePos));
    h5create(file_name, '/global/grid/faces/neighbors', size(model.G.faces.neighbors));
    h5create(file_name, '/global/grid/faces/nodes', size(model.G.faces.nodes));
    h5create(file_name, '/global/grid/faces/areas', size(model.G.faces.areas));
    h5create(file_name, '/global/grid/faces/normals', size(model.G.faces.normals));
    h5create(file_name, '/global/grid/faces/centroids', size(model.G.faces.centroids));
    h5create(file_name, '/global/grid/nodes/num', size(model.G.nodes.num));
    h5create(file_name, '/global/grid/nodes/coords', size(model.G.nodes.coords));
    h5create(file_name, '/global/grid/cartDims', size(model.G.cartDims));
    h5create(file_name, '/global/grid/griddim', size(model.G.griddim));

    % Instantiate negative electrode grid properties
    h5create(file_name, '/local/NegativeElectrode/grid/cells/num', size(model.(ne).G.getNumberOfCells()));
    h5create(file_name, '/local/NegativeElectrode/grid/cells/facePos', size(model.(ne).G.cells.facePos));
    h5create(file_name, '/local/NegativeElectrode/grid/cells/indexMap', size(model.(ne).G.cells.indexMap));
    h5create(file_name, '/local/NegativeElectrode/grid/cells/faces', size(model.(ne).G.cells.faces));
    h5create(file_name, '/local/NegativeElectrode/grid/cells/volumes', size(model.(ne).G.getVolumes()));
    h5create(file_name, '/local/NegativeElectrode/grid/cells/centroids', size(model.(ne).G.cells.centroids));
    h5create(file_name, '/local/NegativeElectrode/grid/faces/num', size(model.(ne).G.getNumberOfFaces()));
    h5create(file_name, '/local/NegativeElectrode/grid/faces/nodePos', size(model.(ne).G.faces.nodePos));
    h5create(file_name, '/local/NegativeElectrode/grid/faces/neighbors', size(model.(ne).G.faces.neighbors));
    h5create(file_name, '/local/NegativeElectrode/grid/faces/nodes', size(model.(ne).G.faces.nodes));
    h5create(file_name, '/local/NegativeElectrode/grid/faces/areas', size(model.(ne).G.faces.areas));
    h5create(file_name, '/local/NegativeElectrode/grid/faces/normals', size(model.(ne).G.faces.normals));
    h5create(file_name, '/local/NegativeElectrode/grid/faces/centroids', size(model.(ne).G.faces.centroids));
    h5create(file_name, '/local/NegativeElectrode/grid/nodes/num', size(model.(ne).G.nodes.num));
    h5create(file_name, '/local/NegativeElectrode/grid/nodes/coords', size(model.(ne).G.nodes.coords));
    h5create(file_name, '/local/NegativeElectrode/grid/cartDims', size(model.(ne).G.cartDims));
    h5create(file_name, '/local/NegativeElectrode/grid/griddim', size(model.(ne).G.griddim));

    % Instantiate positive electrode grid properties
    h5create(file_name, '/local/PositiveElectrode/grid/cells/num', size(model.(pe).G.getNumberOfCells()));
    h5create(file_name, '/local/PositiveElectrode/grid/cells/facePos', size(model.(pe).G.cells.facePos));
    h5create(file_name, '/local/PositiveElectrode/grid/cells/indexMap', size(model.(pe).G.cells.indexMap));
    h5create(file_name, '/local/PositiveElectrode/grid/cells/faces', size(model.(pe).G.cells.faces));
    h5create(file_name, '/local/PositiveElectrode/grid/cells/volumes', size(model.(pe).G.getVolumes()));
    h5create(file_name, '/local/PositiveElectrode/grid/cells/centroids', size(model.(pe).G.cells.centroids));
    h5create(file_name, '/local/PositiveElectrode/grid/faces/num', size(model.(pe).G.getNumberOfFaces()));
    h5create(file_name, '/local/PositiveElectrode/grid/faces/nodePos', size(model.(pe).G.faces.nodePos));
    h5create(file_name, '/local/PositiveElectrode/grid/faces/neighbors', size(model.(pe).G.faces.neighbors));
    h5create(file_name, '/local/PositiveElectrode/grid/faces/nodes', size(model.(pe).G.faces.nodes));
    h5create(file_name, '/local/PositiveElectrode/grid/faces/areas', size(model.(pe).G.faces.areas));
    h5create(file_name, '/local/PositiveElectrode/grid/faces/normals', size(model.(pe).G.faces.normals));
    h5create(file_name, '/local/PositiveElectrode/grid/faces/centroids', size(model.(pe).G.faces.centroids));
    h5create(file_name, '/local/PositiveElectrode/grid/nodes/num', size(model.(pe).G.nodes.num));
    h5create(file_name, '/local/PositiveElectrode/grid/nodes/coords', size(model.(pe).G.nodes.coords));
    h5create(file_name, '/local/PositiveElectrode/grid/cartDims', size(model.(pe).G.cartDims));
    h5create(file_name, '/local/PositiveElectrode/grid/griddim', size(model.(pe).G.griddim));

    % Instantiate electrolyte grid properties
    h5create(file_name, '/local/Electrolyte/grid/cells/num', size(model.(elyte).G.getNumberOfCells()));
    h5create(file_name, '/local/Electrolyte/grid/cells/facePos', size(model.(elyte).G.cells.facePos));
    h5create(file_name, '/local/Electrolyte/grid/cells/indexMap', size(model.(elyte).G.cells.indexMap));
    h5create(file_name, '/local/Electrolyte/grid/cells/faces', size(model.(elyte).G.cells.faces));
    h5create(file_name, '/local/Electrolyte/grid/cells/volumes', size(model.(elyte).G.getVolumes()));
    h5create(file_name, '/local/Electrolyte/grid/cells/centroids', size(model.(elyte).G.cells.centroids));
    h5create(file_name, '/local/Electrolyte/grid/faces/num', size(model.(elyte).G.getNumberOfFaces()));
    h5create(file_name, '/local/Electrolyte/grid/faces/nodePos', size(model.(elyte).G.faces.nodePos));
    h5create(file_name, '/local/Electrolyte/grid/faces/neighbors', size(model.(elyte).G.faces.neighbors));
    h5create(file_name, '/local/Electrolyte/grid/faces/nodes', size(model.(elyte).G.faces.nodes));
    h5create(file_name, '/local/Electrolyte/grid/faces/areas', size(model.(elyte).G.faces.areas));
    h5create(file_name, '/local/Electrolyte/grid/faces/normals', size(model.(elyte).G.faces.normals));
    h5create(file_name, '/local/Electrolyte/grid/faces/centroids', size(model.(elyte).G.faces.centroids));
    h5create(file_name, '/local/Electrolyte/grid/nodes/num', size(model.(elyte).G.nodes.num));
    h5create(file_name, '/local/Electrolyte/grid/nodes/coords', size(model.(elyte).G.nodes.coords));
    h5create(file_name, '/local/Electrolyte/grid/cartDims', size(model.(elyte).G.cartDims));
    h5create(file_name, '/local/Electrolyte/grid/griddim', size(model.(elyte).G.griddim));



    % Write global state variables
    h5write(file_name, '/global/quantities/time', time);
    h5write(file_name, '/global/quantities/voltage', voltage);
    h5write(file_name, '/global/quantities/current', current);

    % Write local state variables
    h5write(file_name, '/local/NegativeElectrode/quantities/surfaceConcentration', ne_surfaceConcentration);
    h5write(file_name, '/local/NegativeElectrode/quantities/concentration', ne_concentration);
    h5write(file_name, '/local/NegativeElectrode/quantities/electricPotential', ne_electricPotential);
    h5write(file_name, '/local/PositiveElectrode/quantities/surfaceConcentration', pe_surfaceConcentration);
    h5write(file_name, '/local/PositiveElectrode/quantities/concentration', pe_concentration);
    h5write(file_name, '/local/PositiveElectrode/quantities/electricPotential', pe_electricPotential);
    h5write(file_name, '/local/Electrolyte/quantities/concentration', elyte_concentration);
    h5write(file_name, '/local/Electrolyte/quantities/electricPotential', elyte_electricPotential);

    % Write specifications
    h5write(file_name, '/global/specification/cRate', CRate);

    % Write global grid properties
    h5write(file_name, '/global/grid/cells/num', model.G.getNumberOfCells());
    h5write(file_name, '/global/grid/cells/facePos', model.G.cells.facePos);
    h5write(file_name, '/global/grid/cells/indexMap', model.G.cells.indexMap);
    h5write(file_name, '/global/grid/cells/faces', model.G.cells.faces);
    h5write(file_name, '/global/grid/cells/volumes', model.G.getVolumes());
    h5write(file_name, '/global/grid/cells/centroids', model.G.cells.centroids);
    h5write(file_name, '/global/grid/faces/num', model.G.getNumberOfFaces());
    h5write(file_name, '/global/grid/faces/nodePos', model.G.faces.nodePos);
    h5write(file_name, '/global/grid/faces/neighbors', model.G.faces.neighbors);
    h5write(file_name, '/global/grid/faces/nodes', model.G.faces.nodes);
    h5write(file_name, '/global/grid/faces/areas', model.G.faces.areas);
    h5write(file_name, '/global/grid/faces/normals', model.G.faces.normals);
    h5write(file_name, '/global/grid/faces/centroids', model.G.faces.centroids);
    h5write(file_name, '/global/grid/nodes/num', model.G.nodes.num);
    h5write(file_name, '/global/grid/nodes/coords', model.G.nodes.coords);
    h5write(file_name, '/global/grid/cartDims', model.G.cartDims);
    h5write(file_name, '/global/grid/griddim', model.G.griddim);

    % Write negative electrode grid properties
    h5write(file_name, '/local/NegativeElectrode/grid/cells/num', model.(ne).G.getNumberOfCells());
    h5write(file_name, '/local/NegativeElectrode/grid/cells/facePos', model.(ne).G.cells.facePos);
    h5write(file_name, '/local/NegativeElectrode/grid/cells/indexMap', model.(ne).G.cells.indexMap);
    h5write(file_name, '/local/NegativeElectrode/grid/cells/faces', model.(ne).G.cells.faces);
    h5write(file_name, '/local/NegativeElectrode/grid/cells/volumes', model.(ne).G.getVolumes());
    h5write(file_name, '/local/NegativeElectrode/grid/cells/centroids', model.(ne).G.cells.centroids);
    h5write(file_name, '/local/NegativeElectrode/grid/faces/num', model.(ne).G.getNumberOfFaces());
    h5write(file_name, '/local/NegativeElectrode/grid/faces/nodePos', model.(ne).G.faces.nodePos);
    h5write(file_name, '/local/NegativeElectrode/grid/faces/neighbors', model.(ne).G.faces.neighbors);
    h5write(file_name, '/local/NegativeElectrode/grid/faces/nodes', model.(ne).G.faces.nodes);
    h5write(file_name, '/local/NegativeElectrode/grid/faces/areas', model.(ne).G.faces.areas);
    h5write(file_name, '/local/NegativeElectrode/grid/faces/normals', model.(ne).G.faces.normals);
    h5write(file_name, '/local/NegativeElectrode/grid/faces/centroids', model.(ne).G.faces.centroids);
    h5write(file_name, '/local/NegativeElectrode/grid/nodes/num', model.(ne).G.nodes.num);
    h5write(file_name, '/local/NegativeElectrode/grid/nodes/coords', model.(ne).G.nodes.coords);
    h5write(file_name, '/local/NegativeElectrode/grid/cartDims', model.(ne).G.cartDims);
    h5write(file_name, '/local/NegativeElectrode/grid/griddim', model.(ne).G.griddim);

    % Write positive electrode grid properties
    h5write(file_name, '/local/PositiveElectrode/grid/cells/num', model.(pe).G.getNumberOfCells());
    h5write(file_name, '/local/PositiveElectrode/grid/cells/facePos', model.(pe).G.cells.facePos);
    h5write(file_name, '/local/PositiveElectrode/grid/cells/indexMap', model.(pe).G.cells.indexMap);
    h5write(file_name, '/local/PositiveElectrode/grid/cells/faces', model.(pe).G.cells.faces);
    h5write(file_name, '/local/PositiveElectrode/grid/cells/volumes', model.(pe).G.getVolumes());
    h5write(file_name, '/local/PositiveElectrode/grid/cells/centroids', model.(pe).G.cells.centroids);
    h5write(file_name, '/local/PositiveElectrode/grid/faces/num', model.(pe).G.getNumberOfFaces());
    h5write(file_name, '/local/PositiveElectrode/grid/faces/nodePos', model.(pe).G.faces.nodePos);
    h5write(file_name, '/local/PositiveElectrode/grid/faces/neighbors', model.(pe).G.faces.neighbors);
    h5write(file_name, '/local/PositiveElectrode/grid/faces/nodes', model.(pe).G.faces.nodes);
    h5write(file_name, '/local/PositiveElectrode/grid/faces/areas', model.(pe).G.faces.areas);
    h5write(file_name, '/local/PositiveElectrode/grid/faces/normals', model.(pe).G.faces.normals);
    h5write(file_name, '/local/PositiveElectrode/grid/faces/centroids', model.(pe).G.faces.centroids);
    h5write(file_name, '/local/PositiveElectrode/grid/nodes/num', model.(pe).G.nodes.num);
    h5write(file_name, '/local/PositiveElectrode/grid/nodes/coords', model.(pe).G.nodes.coords);
    h5write(file_name, '/local/PositiveElectrode/grid/cartDims', model.(pe).G.cartDims);
    h5write(file_name, '/local/PositiveElectrode/grid/griddim', model.(pe).G.griddim);

    % Write electrolyte grid properties
    h5write(file_name, '/local/Electrolyte/grid/cells/num', model.(elyte).G.getNumberOfCells());
    h5write(file_name, '/local/Electrolyte/grid/cells/facePos', model.(elyte).G.cells.facePos);
    h5write(file_name, '/local/Electrolyte/grid/cells/indexMap', model.(elyte).G.cells.indexMap);
    h5write(file_name, '/local/Electrolyte/grid/cells/faces', model.(elyte).G.cells.faces);
    h5write(file_name, '/local/Electrolyte/grid/cells/volumes', model.(elyte).G.getVolumes());
    h5write(file_name, '/local/Electrolyte/grid/cells/centroids', model.(elyte).G.cells.centroids);
    h5write(file_name, '/local/Electrolyte/grid/faces/num', model.(elyte).G.getNumberOfFaces());
    h5write(file_name, '/local/Electrolyte/grid/faces/nodePos', model.(elyte).G.faces.nodePos);
    h5write(file_name, '/local/Electrolyte/grid/faces/neighbors', model.(elyte).G.faces.neighbors);
    h5write(file_name, '/local/Electrolyte/grid/faces/nodes', model.(elyte).G.faces.nodes);
    h5write(file_name, '/local/Electrolyte/grid/faces/areas', model.(elyte).G.faces.areas);
    h5write(file_name, '/local/Electrolyte/grid/faces/normals', model.(elyte).G.faces.normals);
    h5write(file_name, '/local/Electrolyte/grid/faces/centroids', model.(elyte).G.faces.centroids);
    h5write(file_name, '/local/Electrolyte/grid/nodes/num', model.(elyte).G.nodes.num);
    h5write(file_name, '/local/Electrolyte/grid/nodes/coords', model.(elyte).G.nodes.coords);
    h5write(file_name, '/local/Electrolyte/grid/cartDims', model.(elyte).G.cartDims);
    h5write(file_name, '/local/Electrolyte/grid/griddim', model.(elyte).G.griddim);

    H5F.close(fid);

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
