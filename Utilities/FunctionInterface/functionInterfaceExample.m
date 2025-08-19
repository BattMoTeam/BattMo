%% Example of the function interface

%%
% We load a parameter set
jsonstruct = parseBattmoJson('ParameterData/MaterialProperties/LFP/LFP_Xu2015.json');

%%
% In this json structure, we have the description of a function for the OCP. In this case, we have a tabulated function
%
disp(jsonstruct.openCircuitPotential)

%%
% We use the function |setupFunction| to parse the input structure
%
[fn, func] = setupFunction(jsonstruct.openCircuitPotential);

%%
% The function returns a function handler and a |Function| object. The function handler can be directly evaluated.
%

disp(func)

%%
% We have access to the tabulated values

val = [func.dataX, func.dataY];
disp(val)

%%
% The function handler can be used directly as any other matlab function. For example, let us plot this function
%
x = linspace(0, 1, 100);
y = fn(x);

figure
plot(x, y);
xlabel('Stoichiometry');
ylabel('Voltage / V');
title('OCP function LFP Xu 2015')

%%
% The interface is the same for a function defined with an other format. Let us load a function using a string format 

jsonstruct = parseBattmoJson('ParameterData/ParameterSets/Chen2020/chen2020_positive_electrode_interface.json');

%%
% In this json structure, we have the description of a function for the OCP. In this case, we have a tabulated function
%
disp(jsonstruct.openCircuitPotential)
disp(jsonstruct.openCircuitPotential.expressions(1))

%%
% We setup the function and obtain a function object of the class :battmo:`FormulaFunction` and the function handler.
%

[fn, func] = setupFunction(jsonstruct.openCircuitPotential);
disp(func)

%%
% We can use the function handler to evaluate the function

fn(0.5)

%%
% Note that the function handler is a short cut for the generic method |eval| of the parent class |Function|
%

func.eval(0.5)

%%
% As previously, we can plot the function

x = linspace(0, 1, 100);
y = fn(x);

figure
plot(x, y);
xlabel('Stoichiometry');
ylabel('Voltage / V');
title('OCP function NMC Chen 2020')

