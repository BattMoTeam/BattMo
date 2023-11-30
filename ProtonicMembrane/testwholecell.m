clear all

mrstModule add ad-core mrst-gui

clear jsonstruct

filename = 'ProtonicMembrane/gas_supply.json';
jsonstruct.GasSupply = parseBattmoJson(filename);
filename = 'ProtonicMembrane/protonicMembrane.json';
jsonstruct.Cell = parseBattmoJson(filename);

paramobj = ProtonicMembraneCellWithGasSupplyInputParams(jsonstruct);


gen = GasSupplyPEMgridGenerator2D();

gen.nxCell      = 10;
gen.nxGasSupply = 10;
gen.lxCell      = 22*micro*meter;
gen.lxGasSupply = 1*milli*meter;

gen.ny = 10;
gen.ly = 1;

paramobj = gen.updateInputParams(paramobj);

doplot = false;

if doplot
    
    close all

    figure('position', [337, 757, 3068, 557])
    plotGrid(paramobj.G)
    plotGrid(paramobj.Cell.G, 'facecolor', 'red')
    plotGrid(paramobj.GasSupply.G, 'facecolor', 'blue')

    plotGrid(paramobj.GasSupply.G, paramobj.GasSupply.couplingTerms{1}.couplingcells);
    plotGrid(paramobj.GasSupply.G, paramobj.GasSupply.couplingTerms{2}.couplingcells);
    plotGrid(paramobj.Cell.G, paramobj.Cell.couplingTerms{1}.couplingcells(:, 2));
    plotGrid(paramobj.Cell.G, paramobj.Cell.couplingTerms{2}.couplingcells(:, 2));
    plotGrid(paramobj.GasSupply.G, paramobj.GasSupply.couplingTerms{3}.couplingcells );

    return
    
end

model = ProtonicMembraneCellWithGasSupply(paramobj);

initstate = model.setupInitialState();




