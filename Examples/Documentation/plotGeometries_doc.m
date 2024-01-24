doplot.illustration1D = true;
doplot.illustration3D = true;
doplot.jellyroll      = true;
doplot.coincell       = true;

if doplot.illustration1D

    % We fake a 1D model

    jsonstruct_material = parseBattmoJson(fullfile('ParameterData','BatteryCellParameters','LithiumIonBatteryCell','lithium_ion_battery_nmc_graphite.json'));
    jsonstruct_material.include_current_collectors = true;

    paramobj = BatteryInputParams(jsonstruct_material);

    gen = FakeBatteryGeneratorP2D_doc();

    paramobj = gen.updateBatteryInputParams(paramobj);

    model = Battery(paramobj);

    plotBatteryGrid(model);
    
end

if doplot.illustration3D

    % We fake a 1D model

    jsonstruct_material = parseBattmoJson(fullfile('ParameterData','BatteryCellParameters','LithiumIonBatteryCell','lithium_ion_battery_nmc_graphite.json'));
    jsonstruct_material.include_current_collectors = true;

    paramobj = BatteryInputParams(jsonstruct_material);

    gen = BatteryGeneratorP4D();

    paramobj = gen.updateBatteryInputParams(paramobj);

    model = Battery(paramobj);

    plotBatteryGrid(model);
    
end



if doplot.jellyroll
    
    jsonstruct_material = parseBattmoJson(fullfile('ParameterData','BatteryCellParameters','LithiumIonBatteryCell','lithium_ion_battery_nmc_graphite.json'));
    jsonstruct_material.include_current_collectors = true;    
    
    % load json struct for geometry
    jsonstruct_geometry = parseBattmoJson('Examples/JsonDataFiles/4680-geometry.json');
    
    jsonstruct = mergeJsonStructs({jsonstruct_material, jsonstruct_geometry});

    [model, inputparams, jsonstruct, gridGenerator] = setupModelFromJson(jsonstruct);

    figure
    plotGrid(model.G, 'edgealpha', 0.1);
    
end


if doplot.coincell
    
    %% Coin cell

    mrstModule add ad-core mrst-gui mpfa upr

    %% Define some shorthand names for simplicity.
    ne      = 'NegativeElectrode';
    pe      = 'PositiveElectrode';
    co      = 'Coating';
    cc      = 'CurrentCollector';
    elyte   = 'Electrolyte';
    thermal = 'ThermalModel';
    ctrl    = 'Control';
    am      = 'ActiveMaterial';
    sep     = 'Separator';

    %% Setup the properties of Li-ion battery materials and cell design
    jsonstruct = parseBattmoJson(fullfile('ParameterData','BatteryCellParameters','LithiumIonBatteryCell','lithium_ion_battery_nmc_graphite.json'));
    jsonstruct.use_thermal = false;
    jsonstruct.include_current_collectors = true;

    inputparams = BatteryInputParams(jsonstruct);

    %% Setup the geometry and grid for the components
    CRdiameter = 20*milli*meter;
    CRthickness = 1.6*milli*meter;

    compDims = table('rownames', {'NegativeCurrentCollector', ...
                                  'NegativeCoating', ...
                                  'Separator', ...
                                  'PositiveCoating', ...
                                  'PositiveCurrentCollector'});
    numComponents = numel(compDims.Row);

    %% Thickness
    compDims.thickness = zeros(numComponents, 1);

    compDims{'PositiveCoating', 'thickness'} = 67*micro*meter;
    compDims{'Separator'      , 'thickness'} = 20*micro*meter;
    compDims{'NegativeCoating', 'thickness'} = 50*micro*meter;

    currentcollectors = {'PositiveCurrentCollector', 'NegativeCurrentCollector'};
    compDims{currentcollectors, 'thickness'} = 0.5*(CRthickness - sum(compDims.thickness));

    %% Diameters
    compDims.diameter = zeros(numComponents, 1);

    compDims{'PositiveCurrentCollector', 'diameter'} = 1;
    compDims{'PositiveCoating'         , 'diameter'} = 0.8;
    compDims{'Separator'               , 'diameter'} = 0.9;
    compDims{'NegativeCoating'         , 'diameter'} = 0.8;
    compDims{'NegativeCurrentCollector', 'diameter'} = 1;
    compDims.diameter = compDims.diameter * CRdiameter;

    %% Construct grid
    numRadial = 7;
    numAngular = 30;
    hz = min(compDims.thickness);
    compDims.numCellLayers = max(2, round(compDims.thickness / hz));
    compDims{currentcollectors, 'numCellLayers'} = ceil(round(0.25 * compDims{currentcollectors, 'numCellLayers'}));

    disp(compDims);

    params = struct('compDims'  , compDims , ...
                    'numRadial' , numRadial, ...
                    'numAngular', numAngular);

    gen = CoinCellBatteryGenerator();

    % Now, we update the inputparams with the properties of the grid.
    inputparams = gen.updateBatteryInputParams(inputparams, params);

    %%  Initialize the battery model.
    % The battery model is initialized by sending inputparams to the Battery class
    % constructor. see :class:`Battery <Battery.Battery>`.
    model = Battery(inputparams);

    plotBatteryGrid(model)

end

