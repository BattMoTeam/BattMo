
%% Example of variation of one parameter: double layer capacitance
% We offer an example of variations of the EIS spectrum of the P2D model
% when the double layer capacitance changes. 


%%
% First we take the original model developped by Chen et al. in 2020


mrstModule add ad-core mrst-gui mpfa agmg linearsolvers
    
jsonstruct_material = parseBattmoJson(fullfile('ParameterData','ParameterSets','Chen2020','chen2020_lithium_ion_battery.json'));
jsonstruct_geometry = parseBattmoJson(fullfile('Examples', 'JsonDataFiles', 'geometryChen.json'));

jsonstruct = mergeStructs({jsonstruct_material, ...
                               jsonstruct_geometry});

ne = 'NegativeElectrode';
co = 'Coating';
am = 'ActiveMaterial';
itf = 'Interface';

includeDoubleLayer = true;

if includeDoubleLayer

    jsonstruct.(ne).(co).(am).(itf).useDoubleLayerCapacity = true;
    jsonstruct.(ne).(co).(am).(itf).doubleLayerCapacitance = 0.2;

end


[model, inputparams, ~] = setupModelFromJson(jsonstruct);

c_ne = 29.866*mol/litre; % initial concentration at negative electrode
c_pe = 17.038*mol/litre; % initial concentration at positive electrode

initstate = initStateChen2020(model, c_ne, c_pe);

%%
% Then we chose a range of value we evaluate the double layer capacitance
b = logspace(-2, 2, 5);                        

frequences = logspace(-4, 4, 200); 


%%
% Finally we plot the comparative results on one chart.
figure;
hold on;


for i = 1:length(b)
    
    fprintf('compute impedance for doubleLayerCapacitance = %.2e ... ', b(i));

    tic
    jsonstruct.(ne).(co).(am).(itf).doubleLayerCapacitance = b(i);      
    [model, inputparams, ~] = setupModelFromJson(jsonstruct);

    options = [];
    options.stateInitialization.initializationSetup = 'given state';
    options.stateInitialization.computeSteadyState = false;

    extrastructs = [];
    extrastructs.initstate = initstate;
    
    impsolv = ImpedanceSolver(inputparams, options, extrastructs);
    
    Z = impsolv.computeImpedance(frequences);

    fprintf('done in %g s\n', toc);
    
    Z_re = real(Z);
    Z_im = imag(Z);
    curve = sprintf('doubleLayerCapacitance = %.2e', b(i));
    
    plot(Z_re, -Z_im, 'DisplayName', curve);
        
end

legend('show');

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
