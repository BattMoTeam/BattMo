mrstModule add ad-core mrst-gui mpfa agmg linearsolvers
    
jsonstruct_material = parseBattmoJson(fullfile('ParameterData','ParameterSets','Chen2020','chen2020_lithium_ion_battery.json'));
jsonstruct_geometry = parseBattmoJson(fullfile('Examples', 'JsonDataFiles', 'geometryChen.json'));

jsonstruct = mergeJsonStructs({jsonstruct_material, ...
                               jsonstruct_geometry});



includeDoubleLayer = true;

if includeDoubleLayer

    jsonstruct.(ne).(co).(am).(itf).useDoubleLayerCapacity = true;
    jsonstruct.(ne).(co).(am).(itf).doubleLayerCapacitance = 0.2;

end


[model, inputparams, ~] = setupModelFromJson(jsonstruct);

c_ne = 29.866*mol/litre; % initial concentration at negative electrode
c_pe = 17.038*mol/litre; % initial concentration at positive electrode

initstate = initStateChen2020(model, c_ne, c_pe);

b = logspace(-2, 2, 20);                        %change here

frequences = logspace(-2, 4, 30); 

<<<<<<< HEAD
figure;
for i = 1:length(b)
    
    jsonstruct.(ne).(co).(am).(itf).doubleLayerCapacitance = b(i);        %change here
=======

for i = 18:20
>>>>>>> e9e4a96ab3cdbc2f0ce14a0e473826216d412978

    fprintf('compute impedance for doubleLayerCapacitance = %.2e ... ', b(i));  %change here

    tic
    % jsonstruct.(ne).(co).(am).(itf).doubleLayerCapacitance = b(i);        %change here
    % [model, inputparams, ~] = setupModelFromJson(jsonstruct);

    inputparams.(ne).(co).(am).(itf).doubleLayerCapacitance = b(i);
    impsolv = ImpedanceSolver(inputparams, 'initstate', initstate, 'computeSteadyState', false);
    
    Z = impsolv.computeImpedance(frequences);

    fprintf('done in %g s\n', toc);
    
    Z_re = real(Z);
    Z_im = imag(Z);
    curve = sprintf('doubleLayerCapacitance = %.2e', b(i));  %change here
    
    plot(Z_re, -Z_im, 'DisplayName', curve);
        
end

legend('show');