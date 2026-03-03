function jsonstruct = getSimulationJsonInput(jsonstruct)
    
    output = runBatteryJson(jsonstruct, 'runSimulation', false);
    jsonstruct = output.jsonstruct;
    
end

