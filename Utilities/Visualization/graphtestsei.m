%% run stand-alone active material model

% clear the workspace and close open figures
clear
close all

%% Import the required modules from MRST
% load MRST modules
mrstModule add ad-core mrst-gui mpfa

addpath(genpath('/home/xavier/Programs/matlab_bgl/'));

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

paramobj = SingleParticleSEIInputParams(jsonstruct);

paramobj.(an).(sd).N   = 10;
paramobj.(an).(sd).np  = 1;
paramobj.(an).(sei).N  = 10;
paramobj.(an).(sei).np = 1;

xlength = 57e-6; 
G = cartGrid(1, xlength);
G = computeGeometry(G);
paramobj.(an).G = G;

paramobj.(ct).(sd).N   = 10;
paramobj.(ct).(sd).np  = 1;
paramobj.(ct).G = G;

model = SingleParticleSEI(paramobj);

model = model.registerVarAndPropfuncNames();
[g, edgelabels] = setupGraph(model);

dograph = true;
if dograph
    cgt = ComputationalGraphTool(model);
    % cgt.includeNodeNames = 'cSurface';
    cgt.includeNodeNames = 'phiElectrolyte';
    % gg = cgt.setupDescendantGraph();
    % gg = cgt.getComputationalGraph('oneParentOnly', true);    
    [g, edgelabels] = cgt.getComputationalGraph();
    figure
    h = plot(g, 'edgelabel', edgelabels, 'nodefontsize', 10);
end

A = cgt.A;
nodenames = cgt.nodenames;

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
        mn = join(mn, '.');
        if ~isempty(mn)
            mn = mn{1};
            statename = sprintf('state.%s', mn);
        else
            statename = 'state';
        end
        fnname = func2str(fn);
        fnname = regexp(fnname, "\.(.*)", 'tokens');
        fnname = fnname{1}{1};
        fnname = horzcat(mn, {fnname});
        fnname = join(fnname, '.');
        fnname = fnname{1};

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
