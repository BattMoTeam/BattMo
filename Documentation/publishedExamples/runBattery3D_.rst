<?xml version="1.0" encoding="utf-8"?>
<mscript xmlns:mwsh="http://www.mathworks.com/namespace/mcode/v1/syntaxhighlight.dtd">
   <version>9.11</version>
   <release>2021b</release>
   <date>2022-02-16</date>
   <cell>
      <count>1</count>
      <steptitle>Battery 1D model</steptitle>
      <text/>
      <mcode>% load MRST modules
mrstModule add ad-core multimodel mrst-gui mpfa</mcode>
      <mcode-xmlized>
         <mwsh:code xml:space="preserve"><mwsh:comments xml:space="preserve">% load MRST modules</mwsh:comments>
mrstModule <mwsh:strings xml:space="preserve">add</mwsh:strings> <mwsh:strings xml:space="preserve">ad-core</mwsh:strings> <mwsh:strings xml:space="preserve">multimodel</mwsh:strings> <mwsh:strings xml:space="preserve">mrst-gui</mwsh:strings> <mwsh:strings xml:space="preserve">mpfa</mwsh:strings></mwsh:code>
      </mcode-xmlized>
      <mcode-count>1</mcode-count>
      <cellOutputTarget>1</cellOutputTarget>
   </cell>
   <cell>
      <count>2</count>
      <text>
         <p>We create an instance of :class:`BatteryInputParams <a href="Battery.BatteryInputParams">Battery.BatteryInputParams</a>`. This class is used to initiate the battery simulator and it propagates all the parameters through out the submodels.</p>
      </text>
      <mcode>% The input parameters can be given in json format. The json file is read and used to populate the paramobj object.
jsonstruct = parseBatmoJson('JsonDatas/lithiumbattery.json');

paramobj = BatteryInputParams(jsonstruct);

% Some shortcuts used for the sub-models
ne      = 'NegativeElectrode';
pe      = 'PositiveElectrode';
eac     = 'ElectrodeActiveComponent';
cc      = 'CurrentCollector';
elyte   = 'Electrolyte';
thermal = 'ThermalModel';</mcode>
      <mcode-xmlized>
         <mwsh:code xml:space="preserve"><mwsh:comments xml:space="preserve">% The input parameters can be given in json format. The json file is read and used to populate the paramobj object.</mwsh:comments>
jsonstruct = parseBatmoJson(<mwsh:strings xml:space="preserve">'JsonDatas/lithiumbattery.json'</mwsh:strings>);

paramobj = BatteryInputParams(jsonstruct);

<mwsh:comments xml:space="preserve">% Some shortcuts used for the sub-models</mwsh:comments>
ne      = <mwsh:strings xml:space="preserve">'NegativeElectrode'</mwsh:strings>;
pe      = <mwsh:strings xml:space="preserve">'PositiveElectrode'</mwsh:strings>;
eac     = <mwsh:strings xml:space="preserve">'ElectrodeActiveComponent'</mwsh:strings>;
cc      = <mwsh:strings xml:space="preserve">'CurrentCollector'</mwsh:strings>;
elyte   = <mwsh:strings xml:space="preserve">'Electrolyte'</mwsh:strings>;
thermal = <mwsh:strings xml:space="preserve">'ThermalModel'</mwsh:strings>;</mwsh:code>
      </mcode-xmlized>
      <mcode-count>2</mcode-count>
      <cellOutputTarget>2</cellOutputTarget>
   </cell>
   <cell>
      <count>3</count>
      <steptitle>We setup the battery geometry.</steptitle>
      <text>
         <p>Here, we use a 1D model and the class BatteryGenerator1D already contains the discretization parameters</p>
      </text>
      <mcode>gen = BatteryGenerator1D();
gen.fac = 10;
gen = gen.applyResolutionFactors();

% We update pamobj with grid data
paramobj = gen.updateBatteryInputParams(paramobj);

% In this case, we change some of the values of the paramaters that were given in the json file to other values. This is
% done directly on the object paramobj.
paramobj.(ne).(cc).EffectiveElectricalConductivity = 1e5;
paramobj.(pe).(cc).EffectiveElectricalConductivity = 1e5;
paramobj.(thermal).externalTemperature = paramobj.initT;</mcode>
      <mcode-xmlized>
         <mwsh:code xml:space="preserve">gen = BatteryGenerator1D();
gen.fac = 10;
gen = gen.applyResolutionFactors();

<mwsh:comments xml:space="preserve">% We update pamobj with grid data</mwsh:comments>
paramobj = gen.updateBatteryInputParams(paramobj);

<mwsh:comments xml:space="preserve">% In this case, we change some of the values of the paramaters that were given in the json file to other values. This is</mwsh:comments>
<mwsh:comments xml:space="preserve">% done directly on the object paramobj.</mwsh:comments>
paramobj.(ne).(cc).EffectiveElectricalConductivity = 1e5;
paramobj.(pe).(cc).EffectiveElectricalConductivity = 1e5;
paramobj.(thermal).externalTemperature = paramobj.initT;</mwsh:code>
      </mcode-xmlized>
      <mcode-count>3</mcode-count>
      <cellOutputTarget>3</cellOutputTarget>
   </cell>
   <cell>
      <count>4</count>
      <steptitle>The Battery model is initialized by sending paramobj to the Battery class constructor</steptitle>
      <text>
         <p>see :class:`Battery <a href="Battery.Battery">Battery.Battery</a>`</p>
      </text>
      <mcode>model = Battery(paramobj,'use_thermal',true,'use_solid_diffusion',true);
model.AutoDiffBackend= AutoDiffBackend();</mcode>
      <mcode-xmlized>
         <mwsh:code xml:space="preserve">model = Battery(paramobj,<mwsh:strings xml:space="preserve">'use_thermal'</mwsh:strings>,true,<mwsh:strings xml:space="preserve">'use_solid_diffusion'</mwsh:strings>,true);
model.AutoDiffBackend= AutoDiffBackend();</mwsh:code>
      </mcode-xmlized>
      <mcode-count>4</mcode-count>
      <cellOutputTarget>4</cellOutputTarget>
   </cell>
   <cell>
      <count>5</count>
      <steptitle>We compute the cell capacity and chose a discharge rate</steptitle>
      <mcode>C      = computeCellCapacity(model);
CRate  = 1;
inputI = (C/hour)*CRate; % current</mcode>
      <mcode-xmlized>
         <mwsh:code xml:space="preserve">C      = computeCellCapacity(model);
CRate  = 1;
inputI = (C/hour)*CRate; <mwsh:comments xml:space="preserve">% current</mwsh:comments></mwsh:code>
      </mcode-xmlized>
      <mcode-count>5</mcode-count>
      <cellOutputTarget>5</cellOutputTarget>
      <mcodeoutput class="codeoutput">Warning: Adding packing mass in computation of optimal energy 
</mcodeoutput>
   </cell>
   <cell>
      <count>6</count>
      <steptitle>We setup the schedule</steptitle>
      <text>
         <p>We use different time step for the activation phase (small time steps) and the following discharging phase</p>
      </text>
      <mcode>% We start with rampup time steps to go through the activation phase
fac=2;
total = 1.4*hour/CRate;
n=10;
dt0=total*1e-6;
times = getTimeSteps(dt0,n, total,fac);
dt= diff(times);
step = struct('val',diff(times),'control',ones(size(dt)));


% We set up a stopping function. Here, the simulation will stop if the output voltage reach a value smaller than 2. This
% stopping function will not be triggered in this case as we switch to voltage control when E=3.6 (see value of inputE
% below).
pe = 'PositiveElectrode';
cc = 'CurrentCollector';
stopFunc = @(model, state, state_prev) (state.(pe).(cc).E &lt; 2.0);

tup = 0.1; % rampup value for the current function, see rampupSwitchControl
inputE = 3.0; % Value when current control switches to voltage control
srcfunc = @(time, I, E) rampupSwitchControl(time, tup, I, E, inputI, inputE);

% we setup the control by assigning a source and stop function.
control = repmat(struct('src', srcfunc, 'stopFunction', stopFunc), 1, 1);

% This control is used to set up the schedule
schedule = struct('control', control, 'step', step);</mcode>
      <mcode-xmlized>
         <mwsh:code xml:space="preserve"><mwsh:comments xml:space="preserve">% We start with rampup time steps to go through the activation phase</mwsh:comments>
fac=2;
total = 1.4*hour/CRate;
n=10;
dt0=total*1e-6;
times = getTimeSteps(dt0,n, total,fac);
dt= diff(times);
step = struct(<mwsh:strings xml:space="preserve">'val'</mwsh:strings>,diff(times),<mwsh:strings xml:space="preserve">'control'</mwsh:strings>,ones(size(dt)));


<mwsh:comments xml:space="preserve">% We set up a stopping function. Here, the simulation will stop if the output voltage reach a value smaller than 2. This</mwsh:comments>
<mwsh:comments xml:space="preserve">% stopping function will not be triggered in this case as we switch to voltage control when E=3.6 (see value of inputE</mwsh:comments>
<mwsh:comments xml:space="preserve">% below).</mwsh:comments>
pe = <mwsh:strings xml:space="preserve">'PositiveElectrode'</mwsh:strings>;
cc = <mwsh:strings xml:space="preserve">'CurrentCollector'</mwsh:strings>;
stopFunc = @(model, state, state_prev) (state.(pe).(cc).E &lt; 2.0);

tup = 0.1; <mwsh:comments xml:space="preserve">% rampup value for the current function, see rampupSwitchControl</mwsh:comments>
inputE = 3.0; <mwsh:comments xml:space="preserve">% Value when current control switches to voltage control</mwsh:comments>
srcfunc = @(time, I, E) rampupSwitchControl(time, tup, I, E, inputI, inputE);

<mwsh:comments xml:space="preserve">% we setup the control by assigning a source and stop function.</mwsh:comments>
control = repmat(struct(<mwsh:strings xml:space="preserve">'src'</mwsh:strings>, srcfunc, <mwsh:strings xml:space="preserve">'stopFunction'</mwsh:strings>, stopFunc), 1, 1);

<mwsh:comments xml:space="preserve">% This control is used to set up the schedule</mwsh:comments>
schedule = struct(<mwsh:strings xml:space="preserve">'control'</mwsh:strings>, control, <mwsh:strings xml:space="preserve">'step'</mwsh:strings>, step);</mwsh:code>
      </mcode-xmlized>
      <mcode-count>6</mcode-count>
      <cellOutputTarget>6</cellOutputTarget>
   </cell>
   <cell>
      <count>7</count>
      <steptitle>We setup the initial state</steptitle>
      <mcode>initstate = model.setupInitialState();

% Setup nonlinear solver
nls = NonLinearSolver();
% Change default maximum iteration number in nonlinear solver
nls.maxIterations = 10;
% Change default behavior of nonlinear solver, in case of error
nls.errorOnFailure = false;
nls.timeStepSelector=StateChangeTimeStepSelector('TargetProps',{{'PositiveElectrode','CurrentCollector','E'}},'targetChangeAbs',0.03);
% Change default tolerance for nonlinear solver
model.nonlinearTolerance = 1e-3*inputI;
% Set verbosity
model.verbose = true;

% Run simulation

[wellSols, states, report] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls);</mcode>
      <mcode-xmlized>
         <mwsh:code xml:space="preserve">initstate = model.setupInitialState();

<mwsh:comments xml:space="preserve">% Setup nonlinear solver</mwsh:comments>
nls = NonLinearSolver();
<mwsh:comments xml:space="preserve">% Change default maximum iteration number in nonlinear solver</mwsh:comments>
nls.maxIterations = 10;
<mwsh:comments xml:space="preserve">% Change default behavior of nonlinear solver, in case of error</mwsh:comments>
nls.errorOnFailure = false;
nls.timeStepSelector=StateChangeTimeStepSelector(<mwsh:strings xml:space="preserve">'TargetProps'</mwsh:strings>,{{<mwsh:strings xml:space="preserve">'PositiveElectrode'</mwsh:strings>,<mwsh:strings xml:space="preserve">'CurrentCollector'</mwsh:strings>,<mwsh:strings xml:space="preserve">'E'</mwsh:strings>}},<mwsh:strings xml:space="preserve">'targetChangeAbs'</mwsh:strings>,0.03);
<mwsh:comments xml:space="preserve">% Change default tolerance for nonlinear solver</mwsh:comments>
model.nonlinearTolerance = 1e-3*inputI;
<mwsh:comments xml:space="preserve">% Set verbosity</mwsh:comments>
model.verbose = true;

<mwsh:comments xml:space="preserve">% Run simulation</mwsh:comments>

[wellSols, states, report] = simulateScheduleAD(initstate, model, schedule, <mwsh:strings xml:space="preserve">'OutputMinisteps'</mwsh:strings>, true, <mwsh:strings xml:space="preserve">'NonLinearSolver'</mwsh:strings>, nls);</mwsh:code>
      </mcode-xmlized>
      <mcode-count>7</mcode-count>
      <cellOutputTarget>7</cellOutputTarget>
      <mcodeoutput class="codeoutput">Solving timestep 01/18:                                          -&gt; 3 Seconds, 937 Milliseconds
============================================================================================================================================================================================================================================================================================================================================
| It # | elyte_massCons (cell) | elyte_chargeCons (cell) | ne_eac_massCons (cell) | ne_eac_chargeCons (cell) | ne_eac_am_soliddiffeq (sdiff) | pe_eac_massCons (cell) | pe_eac_chargeCons (cell) | pe_eac_am_soliddiffeq (cdiff) | ne_cc_chargeCons (cell) | pe_cc_chargeCons (cell) | energyCons (cell) | EIeq (cell) | controlEq (cntrl) |
============================================================================================================================================================================================================================================================================================================================================
|    1 |*0.00e+00              |*0.00e+00                |*0.00e+00               |*0.00e+00                 |*0.00e+00                      |*0.00e+00               |*0.00e+00                 |*0.00e+00                      |*0.00e+00                |*0.00e+00                |*0.00e+00          |*0.00e+00    | 3.11e+01          |
|    2 |*1.25e-02              |*1.29e-02                |*1.33e-02               |*1.33e-02                 |*2.40e-02                      |*1.09e-02               |*1.09e-02                 | 6.26e-02                      |*1.99e-12                |*5.92e-04                |*1.32e-02          |*5.20e-04    |*0.00e+00          |
|    3 |*8.38e-06              |*1.23e-05                |*1.21e-05               |*1.21e-05                 |*2.25e-05                      |*3.72e-06               |*3.70e-06                 |*2.26e-05                      |*3.30e-13                |*5.92e-04                |*1.68e-04          |*5.20e-04    |*0.00e+00          |
============================================================================================================================================================================================================================================================================================================================================
...</mcodeoutput>
   </cell>
   <cell>
      <count>8</count>
      <steptitle>We process output and recover the output voltage and current from the output states.</steptitle>
      <mcode>ind = cellfun(@(x) not(isempty(x)), states);
states = states(ind);
Enew = cellfun(@(x) x.(pe).(cc).E, states);
Inew = cellfun(@(x) x.(pe).(cc).I, states);
Tmax = cellfun(@(x) max(x.ThermalModel.T), states);
[SOCN,SOCP] =  cellfun(@(x) model.calculateSOC(x), states);
time = cellfun(@(x) x.time, states);</mcode>
      <mcode-xmlized>
         <mwsh:code xml:space="preserve">ind = cellfun(@(x) not(isempty(x)), states);
states = states(ind);
Enew = cellfun(@(x) x.(pe).(cc).E, states);
Inew = cellfun(@(x) x.(pe).(cc).I, states);
Tmax = cellfun(@(x) max(x.ThermalModel.T), states);
[SOCN,SOCP] =  cellfun(@(x) model.calculateSOC(x), states);
time = cellfun(@(x) x.time, states);</mwsh:code>
      </mcode-xmlized>
      <mcode-count>8</mcode-count>
      <cellOutputTarget>8</cellOutputTarget>
   </cell>
   <cell>
      <count>9</count>
      <steptitle>We plot the the output voltage and current</steptitle>
      <mcode>figure
plot((time/hour), Enew, '*-', 'linewidth', 3)
title('Potential (E)')
xlabel('time (hours)')

figure
plot((time/hour), Inew, '*-', 'linewidth', 3)
title('Current (I)')
xlabel('time (hours)')

figure
plot((time/hour), Tmax, '*-', 'linewidth', 3)
title('max(T)')
xlabel('time (hours)')

figure
plot((time/hour), [SOCP,SOCN], '*-', 'linewidth', 3)
title('SOC')
xlabel('time (hours)')
legend('SOC positive','SOC negative')



%{
Copyright 2009-2021 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics &amp; Cybernetics.

This file is part of The Battery Modeling Toolbox BatMo

BatMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BatMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BatMo.  If not, see &lt;http://www.gnu.org/licenses/&gt;.
%}</mcode>
      <mcode-xmlized>
         <mwsh:code xml:space="preserve">figure
plot((time/hour), Enew, <mwsh:strings xml:space="preserve">'*-'</mwsh:strings>, <mwsh:strings xml:space="preserve">'linewidth'</mwsh:strings>, 3)
title(<mwsh:strings xml:space="preserve">'Potential (E)'</mwsh:strings>)
xlabel(<mwsh:strings xml:space="preserve">'time (hours)'</mwsh:strings>)

figure
plot((time/hour), Inew, <mwsh:strings xml:space="preserve">'*-'</mwsh:strings>, <mwsh:strings xml:space="preserve">'linewidth'</mwsh:strings>, 3)
title(<mwsh:strings xml:space="preserve">'Current (I)'</mwsh:strings>)
xlabel(<mwsh:strings xml:space="preserve">'time (hours)'</mwsh:strings>)

figure
plot((time/hour), Tmax, <mwsh:strings xml:space="preserve">'*-'</mwsh:strings>, <mwsh:strings xml:space="preserve">'linewidth'</mwsh:strings>, 3)
title(<mwsh:strings xml:space="preserve">'max(T)'</mwsh:strings>)
xlabel(<mwsh:strings xml:space="preserve">'time (hours)'</mwsh:strings>)

figure
plot((time/hour), [SOCP,SOCN], <mwsh:strings xml:space="preserve">'*-'</mwsh:strings>, <mwsh:strings xml:space="preserve">'linewidth'</mwsh:strings>, 3)
title(<mwsh:strings xml:space="preserve">'SOC'</mwsh:strings>)
xlabel(<mwsh:strings xml:space="preserve">'time (hours)'</mwsh:strings>)
legend(<mwsh:strings xml:space="preserve">'SOC positive'</mwsh:strings>,<mwsh:strings xml:space="preserve">'SOC negative'</mwsh:strings>)



<mwsh:comments xml:space="preserve">%{
</mwsh:comments><mwsh:comments xml:space="preserve">Copyright 2009-2021 SINTEF Industry, Sustainable Energy Technology
</mwsh:comments><mwsh:comments xml:space="preserve">and SINTEF Digital, Mathematics &amp; Cybernetics.
</mwsh:comments><mwsh:comments xml:space="preserve">
</mwsh:comments><mwsh:comments xml:space="preserve">This file is part of The Battery Modeling Toolbox BatMo
</mwsh:comments><mwsh:comments xml:space="preserve">
</mwsh:comments><mwsh:comments xml:space="preserve">BatMo is free software: you can redistribute it and/or modify
</mwsh:comments><mwsh:comments xml:space="preserve">it under the terms of the GNU General Public License as published by
</mwsh:comments><mwsh:comments xml:space="preserve">the Free Software Foundation, either version 3 of the License, or
</mwsh:comments><mwsh:comments xml:space="preserve">(at your option) any later version.
</mwsh:comments><mwsh:comments xml:space="preserve">
</mwsh:comments><mwsh:comments xml:space="preserve">BatMo is distributed in the hope that it will be useful,
</mwsh:comments><mwsh:comments xml:space="preserve">but WITHOUT ANY WARRANTY; without even the implied warranty of
</mwsh:comments><mwsh:comments xml:space="preserve">MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
</mwsh:comments><mwsh:comments xml:space="preserve">GNU General Public License for more details.
</mwsh:comments><mwsh:comments xml:space="preserve">
</mwsh:comments><mwsh:comments xml:space="preserve">You should have received a copy of the GNU General Public License
</mwsh:comments><mwsh:comments xml:space="preserve">along with BatMo.  If not, see &lt;http://www.gnu.org/licenses/&gt;.
</mwsh:comments><mwsh:comments xml:space="preserve">%}</mwsh:comments></mwsh:code>
      </mcode-xmlized>
      <mcode-count>9</mcode-count>
      <cellOutputTarget>9</cellOutputTarget>
      <img src="runBattery1D_01.png"/>
      <img src="runBattery1D_02.png"/>
      <img src="runBattery1D_03.png"/>
      <img src="runBattery1D_04.png"/>
   </cell>
   <originalCode>%% Battery 1D model
% 

% load MRST modules
mrstModule add ad-core multimodel mrst-gui mpfa

%%
% We create an instance of :class:`BatteryInputParams &lt;Battery.BatteryInputParams&gt;`. This class is used to initiate the
% battery simulator and it propagates all the parameters through out the submodels.

% The input parameters can be given in json format. The json file is read and used to populate the paramobj object.
jsonstruct = parseBatmoJson('JsonDatas/lithiumbattery.json');

paramobj = BatteryInputParams(jsonstruct);

% Some shortcuts used for the sub-models
ne      = 'NegativeElectrode';
pe      = 'PositiveElectrode';
eac     = 'ElectrodeActiveComponent';
cc      = 'CurrentCollector';
elyte   = 'Electrolyte';
thermal = 'ThermalModel';

%% We setup the battery geometry.
% Here, we use a 1D model and the class BatteryGenerator1D already contains the discretization parameters
gen = BatteryGenerator1D();
gen.fac = 10;
gen = gen.applyResolutionFactors();

% We update pamobj with grid data
paramobj = gen.updateBatteryInputParams(paramobj);

% In this case, we change some of the values of the paramaters that were given in the json file to other values. This is
% done directly on the object paramobj.
paramobj.(ne).(cc).EffectiveElectricalConductivity = 1e5;
paramobj.(pe).(cc).EffectiveElectricalConductivity = 1e5;
paramobj.(thermal).externalTemperature = paramobj.initT;

%%  The Battery model is initialized by sending paramobj to the Battery class constructor 
% see :class:`Battery &lt;Battery.Battery&gt;`

model = Battery(paramobj,'use_thermal',true,'use_solid_diffusion',true);
model.AutoDiffBackend= AutoDiffBackend();

%% We compute the cell capacity and chose a discharge rate
C      = computeCellCapacity(model);
CRate  = 1; 
inputI = (C/hour)*CRate; % current 

%% We setup the schedule 
% We use different time step for the activation phase (small time steps) and the following discharging phase

% We start with rampup time steps to go through the activation phase 
fac=2;
total = 1.4*hour/CRate;
n=10;
dt0=total*1e-6;
times = getTimeSteps(dt0,n, total,fac);
dt= diff(times);
step = struct('val',diff(times),'control',ones(size(dt)));


% We set up a stopping function. Here, the simulation will stop if the output voltage reach a value smaller than 2. This
% stopping function will not be triggered in this case as we switch to voltage control when E=3.6 (see value of inputE
% below).
pe = 'PositiveElectrode';
cc = 'CurrentCollector';
stopFunc = @(model, state, state_prev) (state.(pe).(cc).E &lt; 2.0); 

tup = 0.1; % rampup value for the current function, see rampupSwitchControl
inputE = 3.0; % Value when current control switches to voltage control
srcfunc = @(time, I, E) rampupSwitchControl(time, tup, I, E, inputI, inputE);

% we setup the control by assigning a source and stop function.
control = repmat(struct('src', srcfunc, 'stopFunction', stopFunc), 1, 1); 

% This control is used to set up the schedule
schedule = struct('control', control, 'step', step); 

%%  We setup the initial state
initstate = model.setupInitialState(); 

% Setup nonlinear solver 
nls = NonLinearSolver(); 
% Change default maximum iteration number in nonlinear solver
nls.maxIterations = 10; 
% Change default behavior of nonlinear solver, in case of error
nls.errorOnFailure = false; 
nls.timeStepSelector=StateChangeTimeStepSelector('TargetProps',{{'PositiveElectrode','CurrentCollector','E'}},'targetChangeAbs',0.03);
% Change default tolerance for nonlinear solver
model.nonlinearTolerance = 1e-3*inputI; 
% Set verbosity
model.verbose = true;

% Run simulation

[wellSols, states, report] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls); 

%%  We process output and recover the output voltage and current from the output states.
ind = cellfun(@(x) not(isempty(x)), states); 
states = states(ind);
Enew = cellfun(@(x) x.(pe).(cc).E, states); 
Inew = cellfun(@(x) x.(pe).(cc).I, states);
Tmax = cellfun(@(x) max(x.ThermalModel.T), states);
[SOCN,SOCP] =  cellfun(@(x) model.calculateSOC(x), states);
time = cellfun(@(x) x.time, states); 

%% We plot the the output voltage and current

figure
plot((time/hour), Enew, '*-', 'linewidth', 3)
title('Potential (E)')
xlabel('time (hours)')

figure
plot((time/hour), Inew, '*-', 'linewidth', 3)
title('Current (I)')
xlabel('time (hours)')

figure
plot((time/hour), Tmax, '*-', 'linewidth', 3)
title('max(T)')
xlabel('time (hours)')

figure
plot((time/hour), [SOCP,SOCN], '*-', 'linewidth', 3)
title('SOC')
xlabel('time (hours)')
legend('SOC positive','SOC negative')



%{
Copyright 2009-2021 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics &amp; Cybernetics.

This file is part of The Battery Modeling Toolbox BatMo

BatMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BatMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BatMo.  If not, see &lt;http://www.gnu.org/licenses/&gt;.
%}
</originalCode>
   <m-file>runBattery1D</m-file>
   <filename>/home/xavier/Matlab/Projects/batmo/Examples/runBattery1D.m</filename>
   <outputdir>/home/xavier/Matlab/Projects/batmo/Documentation/utils/../../Examples/../Documentation/publishedExamples</outputdir>
</mscript>