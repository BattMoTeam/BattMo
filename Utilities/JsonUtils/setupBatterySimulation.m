function simsetup = setupBatterySimulation(jsonstruct)

    output = runBatteryJson(jsonstruct, 'runSimulation', false);
    simsetup = output.simsetup;
    
end

