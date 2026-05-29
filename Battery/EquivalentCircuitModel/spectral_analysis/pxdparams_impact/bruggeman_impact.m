
mrstModule add ad-core mrst-gui mpfa agmg linearsolvers
    
jsonstruct_material = parseBattmoJson(fullfile('ParameterData','ParameterSets','Chen2020','chen2020_lithium_ion_battery.json'));
jsonstruct_geometry = parseBattmoJson(fullfile('Examples', 'JsonDataFiles', 'geometryChen.json'));

jsonstruct = mergeJsonStructs({jsonstruct_material, ...
                               jsonstruct_geometry});



includeDoubleLayer = false;

if includeDoubleLayer

    jsonstruct.(ne).(co).(am).(itf).useDoubleLayerCapacity = true;
    jsonstruct.(ne).(co).(am).(itf).doubleLayerCapacitance = 0.2;

end


[model, inputparams, ~] = setupModelFromJson(jsonstruct);

c_ne = 29.866*mol/litre; % initial concentration at negative electrode
c_pe = 17.038*mol/litre; % initial concentration at positive electrode

initstate = initStateChen2020(model, c_ne, c_pe);
b = linspace(0.5,2,10);

frequences = logspace(-2, 4, 30); 
figure;
hold on;
for i = 1:length(b)
    
    inputparams.NegativeElectrode.Coating.bruggemanCoefficient = b(i);
    inputparams.PositiveElectrode.Coating.bruggemanCoefficient = b(i);
    impsolv = ImpedanceSolver(inputparams, 'initstate', initstate, 'computeSteadyState', false);
    
    Z = impsolv.computeImpedance(frequences);
    Z_re = real(Z);
    Z_im = imag(Z);
    curve = sprintf('Bruggeman coefficient = %.2f', b(i));
    plot(Z_re, -Z_im, 'DisplayName', curve);
end

legend('show');