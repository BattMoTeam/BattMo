%%

clear all
close all

doplot.illustration1D  = false;
doplot.illustration3D  = true;
doplot.jellyroll       = true;
doplot.coincell        = true;
doplot.multilayerpouch = true;
dosave = true;
savedir = fullfile(battmoDir, 'Documentation', 'JOSS', 'figs');

%%
if doplot.illustration1D

    % We fake a 1D model

    jsonstruct_material = parseBattmoJson(fullfile('ParameterData','BatteryCellParameters','LithiumIonBatteryCell','lithium_ion_battery_nmc_graphite.json'));
    jsonstruct_material.include_current_collectors = true;

    inputparams = BatteryInputParams(jsonstruct_material);

    gen = FakeBatteryGeneratorP2D_doc();

    inputparams = gen.updateBatteryInputParams(inputparams);

    model = Battery(inputparams);

    plotBatteryGrid(model);

end

%%
if doplot.illustration3D

    % We fake a 1D model

    jsonstruct_material = parseBattmoJson(fullfile('ParameterData','BatteryCellParameters','LithiumIonBatteryCell','lithium_ion_battery_nmc_graphite.json'));
    jsonstruct_material.include_current_collectors = true;

    inputparams = BatteryInputParams(jsonstruct_material);

    gen = BatteryGeneratorP4D();
    gen.zlength = 10 * gen.zlength;

    inputparams = gen.updateBatteryInputParams(inputparams);

    model = Battery(inputparams);

    plotBatteryGrid(model, 'setstyle', false, 'legendlocation', 'east');

    axis equal tight
    view(3)

    % xlabel('x / mm')
    % ylabel('y / mm')
    % zlabel('z / mm')

    % scaleAxisTicks({'X', 'Y', 'Z'}, 1e3);

    axis off

    if dosave
        exportgraphics(gcf, fullfile(savedir, 'illustration3Dgeometry.png'), 'Resolution', 300);
    end

end

%%

if doplot.jellyroll

    %%
    jsonstruct_material = parseBattmoJson(fullfile('ParameterData','BatteryCellParameters','LithiumIonBatteryCell','lithium_ion_battery_nmc_graphite.json'));
    jsonstruct_material.include_current_collectors = true;

    % load json struct for geometry
    jsonstruct_geometry = parseBattmoJson('Examples/JsonDataFiles/4680-geometry.json');

    jsonstruct_geometry.Geometry.nas = 30;
    jsonstruct_geometry.Geometry.nL = 10;

    jsonstruct = mergeJsonStructs({jsonstruct_material, jsonstruct_geometry});

    model = setupModelFromJson(jsonstruct);

    fig1 = figure;
    plotBatteryGrid(model, 'setstyle', false, 'legend', false, 'fig', fig1);

    axis equal tight off
    camlight left

    if dosave
        exportgraphics(gcf, fullfile(savedir, 'jellyroll_overview.png'), 'Resolution', 300);
    end

    %% Zoom inner

    fig2 = figure;
    plotBatteryGrid(model, 'setstyle', false, 'legend', false, 'fig', fig2);
    axis equal tight off

    cam = SetupCamera(model.grid);
    rInner = jsonstruct_geometry.Geometry.rInner;
    cam.target = [rInner; 0; cam.z];
    cam.viewangle = 1.5;
    cam.phi       = 50;
    cam.theta     = 130;
    cam.do();
    camlight left

    if dosave
        exportgraphics(gcf, fullfile(savedir, 'jellyroll_zoominner.png'), 'Resolution', 300);
    end

    % Zoom outer

    fig3 = figure;
    plotBatteryGrid(model, 'setstyle', false, 'legend', false, 'fig', fig3);
    axis equal tight off

    cam = SetupCamera(model.grid);
    rOuter = jsonstruct_geometry.Geometry.rOuter;
    cam.target = [rOuter; 0; cam.z];
    cam.viewangle = 1.5;
    cam.phi       = 50;
    cam.theta     = 130;
    cam.do();
    camlight left

    if dosave
        exportgraphics(gcf, fullfile(savedir, 'jellyroll_zoomouter.png'), 'Resolution', 300);
    end

end


%%

if doplot.multilayerpouch

    jsonstruct_material = parseBattmoJson(fullfile('ParameterData','BatteryCellParameters','LithiumIonBatteryCell','lithium_ion_battery_nmc_graphite.json'));
    jsonstruct_material.include_current_collectors = true;

    % load json struct for geometry
    jsonstruct_geometry = parseBattmoJson('Examples/JsonDataFiles/geometryMultiLayerPouch.json');

    jsonstruct_geometry.Geometry.nLayers = 30;

    jsonstruct = mergeJsonStructs({jsonstruct_material, jsonstruct_geometry});

    [model, inputparams, jsonstruct, gridGenerator] = setupModelFromJson(jsonstruct);

    plotBatteryGrid(model, 'setstyle', false, 'legend', false);

    axis tight

    scaleAxisTicks({'X', 'Y', 'Z'}, 1e3);
    xlabel('x / mm')
    ylabel('y / mm')
    zlabel('z / mm')

    if dosave
        exportgraphics(gcf, fullfile(savedir, 'multilayerpouch_geometry.png'), 'Resolution', 300);
    end

end

%%

if doplot.coincell

    % Coin cell

    mrstModule add ad-core mrst-gui mpfa upr

    % Define some shorthand names for simplicity.
    ne      = 'NegativeElectrode';
    pe      = 'PositiveElectrode';
    co      = 'Coating';
    cc      = 'CurrentCollector';
    elyte   = 'Electrolyte';
    thermal = 'ThermalModel';
    ctrl    = 'Control';
    am      = 'ActiveMaterial';
    sep     = 'Separator';

    % Setup the properties of Li-ion battery materials and cell design
    jsonstruct = parseBattmoJson(fullfile('ParameterData','BatteryCellParameters','LithiumIonBatteryCell','lithium_ion_battery_nmc_graphite.json'));
    jsonstruct.use_thermal = false;
    jsonstruct.include_current_collectors = true;

    inputparams = BatteryInputParams(jsonstruct);

    % Setup the geometry and grid for the components
    CRdiameter = 20*milli*meter;
    CRthickness = 1.6*milli*meter;

    compDims = table('rownames', {'NegativeCurrentCollector', ...
                                  'NegativeCoating', ...
                                  'Separator', ...
                                  'PositiveCoating', ...
                                  'PositiveCurrentCollector'});
    numComponents = numel(compDims.Row);

    % Thickness
    compDims.thickness = zeros(numComponents, 1);

    compDims{'PositiveCoating', 'thickness'} = 67*micro*meter;
    compDims{'Separator'      , 'thickness'} = 20*micro*meter;
    compDims{'NegativeCoating', 'thickness'} = 50*micro*meter;

    currentcollectors = {'PositiveCurrentCollector', 'NegativeCurrentCollector'};
    compDims{currentcollectors, 'thickness'} = 0.5*(CRthickness - sum(compDims.thickness));

    % Diameters
    compDims.diameter = zeros(numComponents, 1);

    compDims{'PositiveCurrentCollector', 'diameter'} = 1;
    compDims{'PositiveCoating'         , 'diameter'} = 0.8;
    compDims{'Separator'               , 'diameter'} = 0.9;
    compDims{'NegativeCoating'         , 'diameter'} = 0.8;
    compDims{'NegativeCurrentCollector', 'diameter'} = 1;
    compDims.diameter = compDims.diameter * CRdiameter;

    % Construct grid
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

    %  Initialize the battery model.
    % The battery model is initialized by sending inputparams to the Battery class
    % constructor. see :class:`Battery <Battery.Battery>`.
    model = Battery(inputparams);

    plotBatteryGrid(model, 'setstyle', false, 'legend', false);

    axis tight

    xlabel('x / mm')
    ylabel('y / mm')
    zlabel('z / mm')

    scaleAxisTicks({'X', 'Y', 'Z'}, 1e3);

    if dosave
        exportgraphics(gcf, fullfile(savedir, 'coincell_geometry.png'), 'Resolution', 300);
    end

end
