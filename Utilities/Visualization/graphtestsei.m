%% run stand-alone active material model

% clear the workspace and close open figures
clear
close all

%% Import the required modules from MRST
% load MRST modules
mrstModule add ad-core mrst-gui mpfa

% addpath(genpath('/home/xavier/Programs/matlab_bgl/'));

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
    cgf = ComputationalGraphFilter(model);
    % cgf.includeNodeNames = 'cSurface';
    cgf.includeNodeNames = 'phiElectrolyte';
    % gg = cgf.setupDescendantGraph();
    % gg = cgf.setupGraph('oneParentOnly', true);    
    [g, edgelabels] = cgf.setupGraph();
    figure
    h = plot(g, 'edgelabel', edgelabels, 'nodefontsize', 10);
end

return

nn = g.numnodes;
nodenames = g.Nodes.Variables;

[c, ia, ic] = unique(edgelabels.ts);

propnames = c;
propindex = edgelabels.ps(ia);
propindex = cell2mat(propindex);

A = adjacency(g, 'weighted');

doprintspecialvariables = true;

if doprintspecialvariables
    fprintf('Root variables \n');
    nodenames(all(A == 0, 1))
    fprintf('Tail variables \n');
    nodenames(all(A' == 0, 1))
end

p = topological_order(A);

nodenames = nodenames(p);

[lia, locb] = ismember(nodenames, propnames);
propindex = propindex(locb(lia));

funcCallList = {};
for ind = 1 : numel(propindex)
    iprop = propindex(ind);
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

[~, ia, ic] = unique(funcCallList, 'first');

ia = sort(ia);
reducedFuncCallList = funcCallList(ia);

% for ind = 1 : numel(funcCallList)
    % fprintf('%s\n', funcCallList{ind});
% end

fprintf('Function call list\n');
for ind = 1 : numel(reducedFuncCallList)
    fprintf('%s\n', reducedFuncCallList{ind});
end
    
