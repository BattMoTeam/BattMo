:orphan:

.. _runElectrolyser_source:

Source code for runElectrolyser
-------------------------------

.. code:: matlab


  %% Alkaline Membrane Electrolyser 
  
  mrstModule add ad-core matlab_bgl
  
  %% Setup input
  % Setup the physical properties for the electrolyser using json input file
  
  jsonstruct= parseBattmoJson('Electrolyser/Parameters/alkalineElectrolyser.json');
  paramobj = ElectrolyserInputParams(jsonstruct);
  
  %%
  % Setup the grids. We consider a 1D model and the specifications can be read from a json input using
  % :code:`setupElectrolyserGridFromJson`.
  
  jsonstruct= parseBattmoJson('Electrolyser/Parameters/electrolysergeometry1d.json');
  paramobj = setupElectrolyserGridFromJson(paramobj, jsonstruct);
  
  %%
  % We define shortcuts for the different submodels.
  inm = 'IonomerMembrane';
  her = 'HydrogenEvolutionElectrode';
  oer = 'OxygenEvolutionElectrode';
  ptl = 'PorousTransportLayer';
  exr = 'ExchangeReaction';
  ctl = 'CatalystLayer';
  
  %% Setup model
  
  model = Electrolyser(paramobj);
  
  %% Setup the initial condition
  % We use the default initial setup implemented in the model
  
  [model, initstate] = model.setupBcAndInitialState();
  
  %% Setup the schedule with the time discretization
  % We run the simulation over 10 hours, increasing the input current linearly in time.
  
  total = 10*hour;
  
  n   = 100;
  dt  = total/n;
  dts = rampupTimesteps(total, dt, 5);
  
  %%
  % We use the function :code:`rampupControl` to increase the current linearly in time
  
  controlI = -3*ampere/(centi*meter)^2; % if negative, O2 and H2 are produced
  tup      = total; 
  srcfunc  = @(time) rampupControl(time, tup, controlI, 'rampupcase', 'linear');
  control  = struct('src', srcfunc);
  
  step = struct('val', dts, 'control', ones(numel(dts), 1));
  schedule = struct('control', control, 'step', step);
  
  %% Setup the non-linear solver
  % We do only minor modifications here from the standard solver
  
  nls = NonLinearSolver();
  nls.verbose = false;
  nls.errorOnFailure = false;
  
  model.verbose = false;
  
  %% Run the simulation
  
  [~, states, report] = simulateScheduleAD(initstate, model, schedule, 'NonLinearSolver', nls, 'OutputMiniSteps', true);
  
  %% Visualize the results
  %
  % The results contain only the primary variables of the system (the unknwons that descrive the state of the system). We
  % use the method :code:`addVariables` to add all the intermediate quantities that are computed to solve the equations
  % but not stored automatically in the result.
  
  for istate = 1 : numel(states)
      states{istate} = model.addVariables(states{istate});
  end
  
  %%
  % We extract the time, voltage and current values for each time step
  
  time = cellfun(@(state) state.time, states);
  E    = cellfun(@(state) state.(oer).(ptl).E, states);
  I    = cellfun(@(state) state.(oer).(ctl).I, states);
  
  %%
  % We plot the results for the voltage and current
  
  set(0, 'defaultlinelinewidth', 3)
  set(0, 'defaultaxesfontsize', 15)
  
  figure
  subplot(2, 1, 1)
  plot(time/hour, E)
  xlabel('time [hour]');
  ylabel('voltage');
  title('Polarisation curve');
  
  subplot(2, 1, 2)
  plot(time/hour, -I/(1/(centi*meter)^2));
  xlabel('time [hour]');
  ylabel('Current [A/cm^2]');
  title('Input current')
  
  %% pH distribution plot
  %
  % We consider the three domains and plot the pH in each of those. We setup the helper structures to iterate over each
  % domain for the plot.
  
  models = {model.(oer).(ptl), ...
            model.(her).(ptl), ...
            model.(inm)};
  
  fields = {{'OxygenEvolutionElectrode', 'PorousTransportLayer', 'concentrations', 2}  , ...
            {'HydrogenEvolutionElectrode', 'PorousTransportLayer', 'concentrations', 2}, ...
            {'IonomerMembrane', 'cOH'}};
  
  h = figure();
  set(h, 'position', [10, 10, 800, 450]);
  hold on
      
  ntime = numel(time);
  times = linspace(1, ntime, 10);
  cmap  = cmocean('deep', 10);
  
  for ifield = 1 : numel(fields)
  
      fd       = fields{ifield};
      submodel = models{ifield};
  
      x    = submodel.G.cells.centroids;
      
      for itimes = 1 : numel(times);
          
          itime = floor(times(itimes));
          % The method :code:`getProp` is used to recover the value from the state structure
          val   = model.getProp(states{itime}, fd);
          pH    = 14 + log10(val/(mol/litre));
  
          % plot of pH for the current submodel.
          plot(x/(milli*meter), pH, 'color', cmap(itimes, :));
          
      end
  
  end
  
  xlabel('x  /  mm');
  ylabel('pH');
  title('pH distribition in electrolyser')
  
  colormap(cmap)
  hColorbar = colorbar;
  caxis([0 3]);
  hTitle = get(hColorbar, 'Title');
  set(hTitle, 'string', 'J (A/cm^2)');
  
  
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

