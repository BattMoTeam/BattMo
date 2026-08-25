function jsonstruct = getSimulationJsonInput(jsonstruct)
    
    output = runBattery(jsonstruct, 'runSimulation', false);
    jsonstruct = output.jsonstruct;
    
end

