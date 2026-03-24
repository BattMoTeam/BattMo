clear all
close all

printer = @(s) disp(jsonencode(s, 'PrettyPrint', true));

experiment = {
    'Rest for 4000 s';
    'Discharge at 1 mA until 3.0 V';
    'Hold at 3.0 V until 1e-4 A';
    'Charge at 1 A until 4.0 V';
    'Rest for 1 hour';
             };

jsonstruct = convertPyBaMMtoJson(experiment);
printer(jsonstruct);

writeJsonStruct(jsonstruct, 'test.json');

schema = fullfile(battmoDir, 'Utilities', 'JsonSchemas', 'GenericControl.schema.json');
validateJsonStruct(jsonstruct, 'useTmpFile', false, 'schema', schema);
