filename = fullfile(battmoDir()        , ...
                    'Utilities'        , ...
                    'FunctionInterface', ...
                    'Examples'         , ...
                    'classmethodfunctionexample.json');

jsonstruct = parseBattmoJson(filename);

[func1, func1obj] = setupFunction(jsonstruct.func1);
[func2, func2obj] = setupFunction(jsonstruct.func2);
[func3, func3obj] = setupFunction(jsonstruct.func3);
[func4, func4obj] = setupFunction(jsonstruct.func4);

x = 3;
y = 4;

fprintf('Evaluation of func1 : %g\n', func1(x));
fprintf('Evaluation of func2 : %g\n', func2(x, y));
fprintf('Evaluation of func3 : %g\n', func3(x));
fprintf('Evaluation of func4 : %g\n', func4(x));
