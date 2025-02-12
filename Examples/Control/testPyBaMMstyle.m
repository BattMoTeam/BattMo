clear all
close all

printer = @(s) disp(jsonencode(s, 'PrettyPrint', true));

experiment = {
    % 'Discharge at 1C for 0.5 hours';
    % 'Discharge at C/20 for 0.5 hours';
    % 'Charge at 0.5 C for 45 minutes';
    'Discharge at 1 A for 0.5 hours';
    'Charge at 200 mA for 45 minutes';
    %'Discharge at 1W for 0.5 hours';
    %'Charge at 200mW for 45 minutes';
    'Rest for 10 minutes';
    'Hold at 1V for 20 seconds';
    %'Charge at 1 C until 4.1V';
    'Hold at 4.1 V until 50mA';
    % 'Hold at 3V until C/50';
    % 'Discharge at C/3 for 2 hours or until 2.5 V';
             };

for itest = 1:numel(experiment)
    json = convertPyBaMMtoJson(experiment{itest});
    disp(experiment{itest});
    printer(json);
end
