jsonstruct_material = parseBattmoJson(fullfile('ParameterData','ParameterSets','Chen2020','chen2020_lithium_ion_battery.json'));
jsonstruct_geometry = parseBattmoJson(fullfile('Examples', 'JsonDataFiles', 'geometryChen.json'));
jsonstruct = mergeJsonStructs({jsonstruct_material, jsonstruct_geometry});
ne = 'NegativeElectrode';
co = 'Coating';
am = 'ActiveMaterial';
itf = 'Interface';

includeDoubleLayer = true;

if includeDoubleLayer
    
    jsonstruct.(ne).(co).(am).(itf).useDoubleLayerCapacity = true;
    jsonstruct.(ne).(co).(am).(itf).doubleLayerCapacitance = 0.2;
    
end

impexp = ImpedanceExplore(jsonstruct);

% delete(findall(0, 'Type', 'Figure'));
% impexp.parameters{1}{1} = 'PositiveElectrode';
% impexp.setupDefaultParameterLegendNames();
impexp.start();
