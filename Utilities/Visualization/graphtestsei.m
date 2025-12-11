%% run stand-alone active material model

% clear the workspace and close open figures
clear
close all

%% Import the required modules from MRST
% load MRST modules
mrstModule add ad-core mrst-gui mpfa

%% Setup the properties of Li-ion battery materials and cell design
jsonstruct = parseBattmoJson('ParameterData/ParameterSets/Safari2009/fullmodel.json');

% Some shorthands used for the sub-models
an    = 'Anode';
ct    = 'Cathode';
sd    = 'SolidDiffusion';
itf   = 'Interface';
sei   = 'SolidElectrodeInterface';
sr    = 'SideReaction';
elyte = 'Electrolyte';

inputparams = SingleParticleSEIInputParams(jsonstruct);

inputparams.(an).(sd).N   = 10;
inputparams.(an).(sd).np  = 1;
inputparams.(an).(sei).N  = 10;
inputparams.(an).(sei).np = 1;

xlength = 57e-6; 
G = cartGrid(1, xlength);
G = computeGeometry(G);
inputparams.(an).G = G;

inputparams.(ct).(sd).N   = 10;
inputparams.(ct).(sd).np  = 1;
inputparams.(ct).G = G;

model = SingleParticleSEI(inputparams);

model = model.registerVarAndPropfuncNames();
[g, edgelabels] = setupGraph(model);

dograph = true;
if dograph
    cgit = ComputationalGraphInteractiveTool(model);
    % cgit.includeNodeNames = 'cSurface';
    cgit.includeNodeNames = 'phiElectrolyte';
    % gg = cgit.setupDescendantGraph();
    % gg = cgit.getComputationalGraph('oneParentOnly', true);    
    [g, edgelabels] = cgit.getComputationalGraph();
    figure
    h = plot(g, 'edgelabel', edgelabels, 'nodefontsize', 10);
end

A = cgit.A;
nodenames = cgit.nodenames;

doprintspecialvariables = true;

if doprintspecialvariables
    fprintf('Root variables \n');
    nodenames(all(A == 0, 1))
    fprintf('Tail variables \n');
    nodenames(all(A' == 0, 1))
end

p = topological_order(A);

funcCallList = {};
for ind = 1 : numel(p)
    iprop = full(A(:, p(ind)));
    iprop = unique(iprop(iprop>0));
    if ~isempty(iprop)
        assert(numel(iprop) == 1, 'problem');
        propfunction = model.propertyFunctionList{iprop};
        fn = propfunction.fn;
        mn = propfunction.modelnamespace;
        mn = strjoin(mn, '.');
        if ~isempty(mn)
            statename = sprintf('state.%s', mn);
        else
            statename = 'state';
        end
        fnname = func2str(fn);
        fnname = regexp(fnname, "\.(.*)", 'tokens');
        fnname = fnname{1}{1};
        fnname = horzcat(mn, {fnname});
        fnname = strjoin(fnname, '.');

        funcCallList{end + 1} = sprintf('%s = model.%s(%s);', statename, fnname, statename);
    end
end

[~, ia, ic] = unique(funcCallList, 'first');

ia = sort(ia);
reducedFuncCallList = funcCallList(ia);

fprintf('Function call list\n');
for ind = 1 : numel(reducedFuncCallList)
    fprintf('%s\n', reducedFuncCallList{ind});
end



%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
