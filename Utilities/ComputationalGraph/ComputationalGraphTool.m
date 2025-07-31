classdef ComputationalGraphTool
%
% A model is essentially characterised by its computational graph. The functions to setup the graph of a model are given in the BaseModel class.
%
% The ComputationalGraphTool is used to store the computational graph and provides utility functions to explore it in an iteractive way.
%
% From a model, you can get the computatonial graph by writing in the terminal
%
% cgt = model.cgt
%
% Then write
%
% cgt.help

% to get an overview of the different commands that are available.
% 
    properties (SetAccess = private)

        model

        functionDocs % cell array. Each cell contains a struct describing the function documentation, with fields
                     % - name      : function name
                     % - docstring : string describing the function

    end

    properties
        
        adjencyMatrix    % Adjency matrix for the computational graph
                         % - column index      : output variable index (as in cgt.varNameList)
                         % - row index         : input variable (as in cgt.varNameList)
                         % - coefficient value : property function index as in model.propertyFunctionList
        varNameList      % List of variable names (object of class VarName).
                         % They have been processed model.varNameList so that there exist separate entries for each index in the case of variable names with indices.
        nodenames        % The node names in the graph are unique string representations of the variable names in varNameList.
        staticprops      % list of the static properties (property functions that take no input). It is computed in setupComputationalGraph.
                         % Each element is a struct with the  the fields
                         % - nodename (element in nodenames)
                         % - propind (index of the property function in the list model.propertyFunctionList)
                         % - varnameind (index of the variable name in varNameList)
        extraVarNameInds % index relative with respect to varNameList
        isok             % The graph is a proper directed acyclic graph that can be used for computation.

        %% model graph (for visualization only, i.e. not used in assembly)
        modelnames         % List of strings given the model names (contructed using getSubModelNames from BaseModel)
        modelAdjencyMatrix % model adjency matrix

    end

    methods

        function cgt = ComputationalGraphTool(model)
            if isempty(model.propertyFunctionList)
                model = model.registerVarAndPropfuncNames();
            end
            cgt.model = model;
            [A, staticprops, varNameList, nodenames] = setupGraph(model);

            % In adjacency matrix A,
            % - column index      : output variable index (as in cgt.varNameList)
            % - row index         : input variable (as in cgt.varNameList)
            % - coefficient value : property function index as in model.propertyFunctionList

            cgt.adjencyMatrix = A;
            cgt.varNameList   = varNameList;
            cgt.nodenames     = nodenames;
            cgt.staticprops   = staticprops;
            cgt.isok        = false;
            if size(A, 1) == numel(varNameList)
                cgt = cgt.setupComputationalGraph();
            else
                fprintf('\nThe graph could not be ordered properly due a mismatch in the variable declarations\nFix that before using the graph in computations\n');
            end

        end

        function cgt = setupComputationalGraph(cgt)

            A            = cgt.adjencyMatrix;
            varNameList  = cgt.varNameList;
            nodenames    = cgt.nodenames;
            staticprops  = cgt.staticprops;

            p = compatible_topological_order(A);

            if isempty(p)
                fprintf('The graph contains cycles. It implies that some variables cannot be evaluated.\n')
                cgt.isok = false;
                return
            end

            varNameList = varNameList(p);
            A           = A(p, p);
            nodenames   = nodenames(p);

            for istat = 1 : numel(staticprops)
                nodename = staticprops{istat}.nodename;
                [isok, varnameind] = ismember(nodename, nodenames);
                if ~isok
                    fprintf('The static variable %s has not been found in list\n', nodename);
                    cgt.isok = false;
                    return
                end

                staticprops{istat}.varnameind = varnameind;
            end

            cgt.varNameList   = varNameList;
            cgt.adjencyMatrix = A;
            cgt.nodenames     = nodenames;
            cgt.staticprops   = staticprops;

            extraVarNames = cgt.model.extraVarNameList;
            extraVarNameInds = [];
            for ivar = 1 : numel(extraVarNames);
                extraVarName = extraVarNames{ivar};
                extraVarName_s = extraVarName.resolveIndex();
                for iivar = 1 : numel(extraVarName_s)
                    extraVarNameInds(end + 1) = cgt.getVarNameIndex(extraVarName_s{iivar});
                end
            end
            cgt.extraVarNameInds = extraVarNameInds;

            cgt.isok = true;
        end

        function [varnames, varnameinds, propfuncinds, levels] = getDependencyList(cgt, varname, direction)
        % dependency direction can be
        % - 'upwards' for upwards in the graph
        % - 'downwards' for downwards in the graph
            
            A = cgt.adjencyMatrix;

            % for cell-valued variable pick-up a valid index
            if (varname.dim > 1) && (ischar(varname.index))
                varname.index = 1;
            elseif isnumeric(varname.index) && (numel(varname.index) > 1)
                varname.index = varname.index(1);
            end

            varnameind = cgt.getVarNameIndex(varname);

            switch direction
              case 'upwards'
                [depvarnameinds, propfuncinds, propvarnameinds, propdeplevels, rootinds, rootdeplevels] = getDependencyVarNameInds(varnameind, A);
              case 'downwards'
                [depvarnameinds, propfuncinds, propvarnameinds, propdeplevels, rootinds, rootdeplevels] = getDependencyVarNameInds(varnameind, A');
              otherwise
                error('direction not recognized')
            end
            
            varnameinds = depvarnameinds;
            levels      = [rootdeplevels; propdeplevels];
            varnames    = cgt.varNameList(varnameinds);

        end

        function printChildDependencyList(cgt, varname)
            
            cgt.printDependencyList(varname, 'downwards');
            
        end
        
        function printParentDependencyList(cgt, varname)
            
            cgt.printDependencyList(varname, 'upwards');

        end

        function varnameind = findVarName(cgt, varname)

            nodenames = cgt.nodenames;

            if isa(varname, 'VarName')
                nodename = ComputationalGraphTool.getNodeName(varname);
                varnameind = findVarName(nodename);
            elseif isa(varname, 'cell')
                varname = VarName(varname(1 : end - 1), varname{end});
                varnameind = findVarName(varname);
            elseif isa(varname, 'char')
                varnameinds = regexpSelect(cgt.nodenames, varname);
                if numel(varnameinds) > 1
                    fprintf('Several variables found:\n\n')
                    for ivar = 1 : numel(varnameinds)
                        fprintf('%d : %s\n', ivar, nodenames{varnameinds(ivar)});
                    end
                    ivar = input('Pick one:');
                    varnameind = varnameinds(ivar);
                else
                    varnameind = varnameinds;
                end
            else
                error('input varname not recognized');
            end

        end
        
        function printDependencyList(cgt, varname, direction)
        % input varname is either
        %  - a VarName instance
        %  - a cell which then uses shortcuts for VarName (see implementation below)
        %  - a string giving a regexp. It will be used to select varnames by the string name
        %
        % dependency direction can be
        % - 'upwards' for upwards in the graph
        % - 'downwards' for downwards in the graph
        %
        % Prints the dependency list of the variable given by varname
            

            varnameind = cgt.findVarName(varname);
            varname = cgt.varNameList{varnameind};
            
            [varnames, varnameinds, propfuncinds, distance] = cgt.getDependencyList(varname, direction);

            for ivar = 1 : numel(varnameinds)
                varnameind = varnameinds(ivar);
                fprintf('%s (%d)\n', cgt.nodenames{varnameind}, distance(ivar));
            end

        end

        function [propfuncs, propfuncinds, varnameinds] = getPropFunctionList(cgt, propfunc)

        % Get the list of property functions and the corresponding indices (with respect to
        % cgt.model.propertyFunctionList) that should be call so that the variable of propfunc get updated.

            A           = cgt.adjencyMatrix;
            staticprops = cgt.staticprops;

            varname = propfunc.varname;

            % for cell-valued variable pick-up a valid index
            if (varname.dim > 1) && (ischar(varname.index))
                varname.index = 1;
            elseif isnumeric(varname.index) && (numel(varname.index) > 1)
                varname.index = varname.index(1);
            end

            varnameind = cgt.getVarNameIndex(varname);

            [~, propfuncinds, propvarnameinds, ~, rootinds, ~] = getDependencyVarNameInds(varnameind, A);

            [rootinds, rootpropinds] = cgt.findStaticVarNameInds(rootinds);

            if ~isempty(rootinds)
                propfuncinds = [rootpropinds; propfuncinds];
                varnameinds  = [rootinds; propvarnameinds];
            else
                varnameinds = propvarnameinds;
            end

            propfuncs = cgt.model.propertyFunctionList(propfuncinds);

        end

        function [staticinds, staticpropinds] = findStaticVarNameInds(cgt, varnameinds)

        % Given a set list of index of varnameinds (index with respect to cgt.varNameList), extract and return the indices that
        % correspond to static variables.

            staticprops = cgt.staticprops;

            allstaticinds     = cellfun(@(staticprop) staticprop.varnameind, staticprops)';
            allstaticpropinds = cellfun(@(staticprop) staticprop.propind, staticprops)';

            [isok, ia] = ismember(varnameinds, allstaticinds);

            if any(isok)
                staticpropinds = allstaticpropinds(ia(isok), 1);
                staticinds     = allstaticinds(ia(isok), 1);
            else
                staticpropinds = {};
                staticinds = {};
            end

        end

        function funcCallList = getPropFunctionCallList(cgt, propfunc)
        % input propfunc is either
        % - an instance of PropFunction
        % - a valid input for findPropFunction, that is, either
        %     - a VarName instance
        %     - a cell which then uses shortcuts for VarName (see implementation in findPropFunction)
        %     - a string giving a regexp. It will be used to select varnames by the string name
        %   In this case, findPropFunction is first run to obtain a list of propfunctions
        % setup the list of function call (as list of cell of strings) that will be run to update the property function propfunc.

            if isa(propfunc, 'PropFunction')
                propfuncs = cgt.getPropFunctionList(propfunc);
            elseif isa(propfunc, 'VarName')
                % We set up a propfunction with only the varname
                propfunc = PropFunction(propfunc, [], [], []);
                propfuncs = cgt.getPropFunctionList(propfunc);
            else

                varname = propfunc;
                selectedpropfuncs = cgt.findPropFunction(propfunc);
                if isa(selectedpropfuncs, 'PropFunction')
                    selectedpropfuncs = {selectedpropfuncs};
                end
                propfuncs   = {};
                for isel = 1 : numel(selectedpropfuncs)
                    deppropfuncs = cgt.getPropFunctionList(selectedpropfuncs{isel});
                    propfuncs   = horzcat(propfuncs, deppropfuncs);
                end
            end

            funcCallList = {};

            for iprop = 1 : numel(propfuncs)

                propfunc = propfuncs{iprop};
                fncallstr = propfunc.functionCallSetupFn(propfunc);
                if ~isempty(fncallstr)
                    funcCallList{end + 1} = fncallstr;
                end

            end

        end

        function printPropFunctionCallList(cgt, propfunc, varargin)
        % input propfunc is either
        % - an instance of PropFunction
        % - a valid input for findPropFunction, that is, either
        %     - a VarName instance
        %     - a cell which then uses shortcuts for VarName (see implementation below)
        %     - a string giving a regexp. It will be used to select varnames by the string name
        %   In this case, findPropFunction is first run to obtain the property function. It should return a unique element, otherwise we get a warning with a list of returned elements
        % Print the list of function call (as string) that will update the property function propfunc.


            opt = struct('fullSignature', false);
            opt = merge_options(opt, varargin{:});

            if ~isa(propfunc, 'PropFunction')

                varname = propfunc;

                propfunc = cgt.findPropFunction(varname);

                if isempty(propfunc)
                    fprintf('No property matching regexp has been found\n');
                    return
                end

                if numel(propfunc) > 1
                    fprintf('Several property functions are matching\n\n');
                    cgt.printPropFunction(varname);
                    return
                end

            end

            nodenames = cgt.nodenames;

            [propfuncs, ~, varnameinds] = cgt.getPropFunctionList(propfunc);

            funcCallList = {};

            strls = arrayfun(@(varnameind) strlength(nodenames{varnameind}), varnameinds);
            strl = max(strls);

            for iprop = 1 : numel(propfuncs)

                propfunc = propfuncs{iprop};
                varnameind = varnameinds(iprop);

                if opt.fullSignature
                    fncallstr = propfunc.functionCallSetupFn(propfunc);
                    varstr = sprintf('state.%s', nodenames{varnameind});
                    varstr = sprintf('%-*s', strl + 6, varstr);
                    fprintf('%s <- %s\n', varstr, fncallstr);
                else
                    inputvarnames = propfunc.inputvarnames;
                    strs = {};
                    for ivar = 1 : numel(inputvarnames)
                        inputvarname = inputvarnames{ivar};
                        strs = horzcat(strs, cgt.getNodeName(inputvarname));
                    end
                    strs = strjoin(strs, ', ');
                    varstr = sprintf('%-*s', strl + 6, nodenames{varnameind});
                    fprintf('%s [%s]\n', varstr, strs);
                end

            end

        end

        function nodestaticnames = getStaticVarNames(cgt)

        % add the variables that are updated with no input arguments
            propfunctions = cgt.model.propertyFunctionList;
            varnames = {};
            for iprop = 1 : numel(propfunctions)
                propfunction = propfunctions{iprop};
                if isempty(propfunction.inputvarnames)
                    varnames{end + 1} = propfunction.varname;
                end
            end

            nodestaticnames = cgt.getNodeNames(varnames);

        end

        function [propfuncs, propfuncinds] = findPropFunction(cgt, varname)
        % The input varnames can be either:
        % - a VarName instance
        % - a cell which then uses shortcuts for VarName (see implementation below)
        % - a string giving a regexp. It will be used to select varnames by the string name
        %
        % The function returns a list of PropFunction that updates the variable given by varname.

            nodenames = cgt.nodenames;
            A         = cgt.adjencyMatrix;
            model     = cgt.model;

            if isa(varname, 'VarName')
                varnameinds = cgt.getVarNameIndex(varname);
            elseif isa(varname, 'cell')
                varname = VarName(varname(1 : end - 1), varname{end});
                [propfuncs, propfuncinds] = cgt.findPropFunction(varname);
                return
            elseif isa(varname, 'char')
                selectednodenames = varname;
                varnameinds = regexpSelect(cgt.nodenames, varname);
            else
                error('input type not recognized');
            end

            propfuncinds = max(A(:, varnameinds), [], 1)';
            staticinds = varnameinds(propfuncinds == 0);
            propfuncinds = propfuncinds(propfuncinds > 0); % remove the zero elements

            [staticinds, staticpropinds] = cgt.findStaticVarNameInds(staticinds);

            if ~isempty(staticinds)
                propfuncinds = [staticpropinds; propfuncinds];
            end

            propfuncinds = unique(propfuncinds);
            
            propfuncs = model.propertyFunctionList(propfuncinds);

            if numel(propfuncs) == 1
                propfuncs = propfuncs{1};
            end

        end


        function openPropFunction(cgt, nodename)
        % Open in editor the place where the variable that matches the regexp nodename is updated. The regexp nodename
        % should return a unique match

            [propfunc, propfuncinds] = cgt.findPropFunction(nodename);

            if isempty(propfunc)
                fprintf('No property matching regexp has been found\n');
                return
            end

            if numel(propfunc) > 1
                fprintf('Several property functions are matching\n\n');
                for iprop = 1 : numel(propfunc)
                    varname = propfunc{iprop}.varname;
                    nodenames = ComputationalGraphTool.getNodeName(varname);
                    for inode = 1 : numel(nodenames)
                        fprintf('%d : %s\n', iprop, nodenames{inode});
                    end
                end
                iprop = input('Pick one:');
                propfunc = propfunc{iprop};
            end

            % Initialise variable state with empty structure. In this way calling the property function will raise an error that
            % we will catch (see call passed in eval below)
            mn = propfunc.modelnamespace;
            
            state         = ComputationalGraphTool.setupState([], mn);
            state0        = state; % may be needed in case of accumulation term
            dt            = 0;     % may be needed in case of accumulation term
            drivingForces = [];    % may be needed in case of update of driving force term
            model         = cgt.model; % needed in function call passed in eval
            
            fncallstr = propfunc.functionCallSetupFn(propfunc);

            try
                % fn = @(model,state)ProtonicMembraneGasSupply.updateDensity(model,state);
                % state.GasSupplyBc = fn(model.GasSupplyBc, state.GasSupplyBc);
                eval(fncallstr);
            catch ME
                fprintf('%s\n', fncallstr);
                stack = ME.stack;
                file = stack.file;
                lineNum = stack.line;
                editor = settings().matlab.editor.OtherEditor.ActiveValue;
                if contains(editor, 'emacs')
                    cmd = sprintf('%s +%d %s &', editor, lineNum, file);
                    system(cmd);
                else
                    matlab.desktop.editor.openAndGoToLine(file, lineNum);
                end
            end

        end

        function printVarNames(cgt, nodename)
        % Print the variable(s) that match the regexp nodename

            nodenames = cgt.nodenames;

            if nargin < 2
                nodename = '.';
            end
            inds = regexpSelect(cgt.nodenames, nodename);
            nodenames = nodenames(inds);
            for ind = 1 : numel(nodenames)
                fprintf('%s\n', nodenames{ind});
            end

        end


        function printPropFunction(cgt, nodename)
        % Print property function including output variable name, function name and input variable names for the variable(s)
        % that match the regexp nodename

            propfuncs = cgt.findPropFunction(nodename);

            if numel(propfuncs) == 1
                propfuncs = {propfuncs};
            end

            for iprop = 1 : numel(propfuncs)

                propfunc      = propfuncs{iprop};
                varname       = propfunc.varname;
                varname_s     = varname.resolveIndex();
                outputvarstrs = {};

                fncallstr = propfunc.getFunctionCallString();

                if ~isempty(fncallstr)
                    for ind = 1 : numel(varname_s)
                        fullname = varname_s{ind}.getIndexedFieldname();
                        outputvarstrs{end + 1} = sprintf('state.%s', fullname);
                    end
                    outputvarstr = strjoin(outputvarstrs, {', '});
                    fprintf('%s <-', outputvarstr);
                    fprintf(' (%s) ', fncallstr);
                    inputvarnames = propfunc.inputvarnames;
                    if isempty(inputvarnames)
                        fprintf('[no state field is used]\n');
                    else
                        fprintf(' <- ', fncallstr);
                        inputvarstrs = {};
                        for ind = 1 : numel(inputvarnames)
                            varname = inputvarnames{ind};
                            inputvarstrs{end + 1} = sprintf('state.%s', varname.getFieldname());
                        end
                        inputvarstr = strjoin(inputvarstrs, {', '});
                        fprintf('[%s]\n', inputvarstr);
                    end
                end
            end

        end

        function varnameind = getVarNameIndex(cgt, varname)
        % Given a varname, find its index in cgt.varNameList
            for varnameind = 1 : numel(cgt.varNameList)
                if varname.compareVarName(cgt.varNameList{varnameind})
                    return
                end
            end

            fprintf('Variable not found\n');
            varnameind = [];


        end

        function cgt = setupModelGraph(cgt, varargin)


            function g = addModel(g, model, modelname, namespace, parentnodename)

                namecomponents = horzcat(namespace, {modelname});
                nodename = strjoin(namecomponents, '.');

                g = addnode(g, nodename);

                if ~isempty(parentnodename)
                    g = addedge(g, parentnodename, nodename);
                end
                submodelnames = model.getSubModelNames();
                namespace = horzcat(namespace, modelname);

                for isubmodel = 1 : numel(submodelnames)

                    submodelname = submodelnames{isubmodel};
                    submodel = model.(submodelname);

                    g = addModel(g, submodel, submodelname, namespace, nodename);

                end

            end

            g = digraph();
            g = addModel(g, cgt.model, 'model', {}, []);

            cgt.modelnames = g.Nodes.Variables;
            cgt.modelAdjencyMatrix = adjacency(g, 'weighted');

        end

        function cgt2 = printModelNames(cgt, modelname)

            if isempty(cgt.modelnames)
                cgt = cgt.setupModelGraph();
                if nargout > 0
                    cgt2 = cgt;
                end
            end

            if nargin < 2
                modelname = '.';
            end

            modelnames = cgt.modelnames;
            inds = regexpSelect(modelnames, modelname);
            modelnames = modelnames(inds);

            for imodel = 1 : numel(modelnames)
                fprintf('%s\n', modelnames{imodel});
            end

        end


        function printRootVariables(cgt, nodename)
        % Print the root variables in computational graph

            A = cgt.adjencyMatrix;
            nodenames = cgt.nodenames;

            %% print root variables after removing the variables that were declared as static in the model

            rootnames = nodenames(all(A == 0, 1));
            staticnames = cgt.getStaticVarNames();

            ind = ismember(rootnames, staticnames);
            rootinds = find(~ind);

            if nargin > 1
                inds = regexpSelect(cgt.nodenames, nodename);
                rootinds = intersect(inds, rootinds);
            end
            
            cgt.printHeader('Root Variables', numel(rootinds));

            for irind = 1 : numel(rootinds)
                fprintf('%s\n', rootnames{rootinds(irind)});
            end

            if nnz(ind)
                staticinds = find(ind);
                if nargin > 1
                    inds = regexpSelect(cgt.nodenames, nodename);
                    staticinds = intersect(inds, staticinds);
                end
                cgt.printHeader('Static Variables', numel(staticinds));
                for isind = 1 : numel(staticinds)
                    fprintf('%s\n', rootnames{staticinds(isind)});
                end
            end

        end

        function printTailVariables(cgt, nodename)
        % Print the tail variables of the computational graph

            A                = cgt.adjencyMatrix;
            nodenames        = cgt.nodenames;
            extravarnameinds = cgt.extraVarNameInds;

            tailinds = find(all(A' == 0, 1));
            isadded  = ismember(tailinds, extravarnameinds);
            tailinds = tailinds(~isadded);

            if nargin > 1
                inds = regexpSelect(cgt.nodenames, nodename);
                tailinds = intersect(inds, tailinds);
            end

            cgt.printHeader('Tail Variables', numel(tailinds));
            
            for itail = 1 : numel(tailinds)
                fprintf('%s\n', nodenames{tailinds(itail)});
            end

            if nargin == 1
                cgt.printHeader('Extra Variables (do not enter in residual evaluation)', numel(extravarnameinds));
                for ievar = 1 : numel(extravarnameinds)
                    fprintf('%s\n', nodenames{extravarnameinds(ievar)});
                end
            end

        end

        function printDetachedVariables(cgt)
        % Print the "detached" variables, which are the variables that are not connected to the graph. This is specially useful
        % in debugging because such variables should be eliminated from the final graph.

            A = cgt.adjencyMatrix;
            nodenames = cgt.nodenames;

            fprintf('Detached variables \n');
            ind1 = all(A == 0, 1);
            ind2 = all(A' == 0, 1);
            detached_nodenames = nodenames(ind1&ind2);

            for idnodename = 1 : numel(detached_nodenames)
                fprintf('%s\n', detached_nodenames(idnodename));
            end
        end

        function primvarnames = getPrimaryVariableNames(cgt)
        % Return the primary variables, which are defined as the root variables and not declared or recognized as static.

            A = cgt.adjencyMatrix;
            nodenames = cgt.nodenames;

            ind = find(all(A == 0, 1));

            rootnames = nodenames(ind);
            staticnames = cgt.getStaticVarNames();
            indstatic = ismember(rootnames, staticnames);
            ind = ind(~indstatic);
            primVarNameList = cgt.varNameList(ind);

            primvarnames = cellfun(@(varname) varname.getPropName(), primVarNameList, 'uniformoutput', false);

        end

        function eqvarnames = getEquationVariableNames(cgt)
        % Return the equation variables, which are defined as the tail variables and not declared as extravarnames.

            A                = cgt.adjencyMatrix;
            nodenames        = cgt.nodenames;
            extravarnameinds = cgt.extraVarNameInds;

            eqinds = find(all(A' == 0, 1));
            isadded  = ismember(eqinds, extravarnameinds);
            eqinds = eqinds(~isadded);

            equationVarNameList = cgt.varNameList(eqinds);

            eqvarnames = cellfun(@(varname) varname.getPropName(), equationVarNameList, 'uniformoutput', false);

        end


        function funcCallList = getOrderedFunctionCallList(cgt, varargin)
        % Return a list of strings that corresponds to all the function calls ordered in the right order.

            opt = struct('removeExtraVariables', true);
            opt = merge_options(opt, varargin{:});

            A                = cgt.adjencyMatrix;
            staticprops      = cgt.staticprops;
            extraVarNameInds = cgt.extraVarNameInds;

            funcCallList = {};

            if ~isempty(staticprops)

                for istatic = 1 : numel(staticprops)

                    iprop = staticprops{istatic}.propind;
                    propfunction = cgt.model.propertyFunctionList{iprop};
                    fncallstr = propfunction.getFunctionCallString();
                    if ~isempty(fncallstr)
                        funcCallList{end + 1} = fncallstr;
                    end

                end

            end

            if opt.removeExtraVariables
                nA = size(A, 2);
                ind = true(nA, 1);
                ind(extraVarNameInds) = false;
                A = A(ind, ind);
            end

            for ind = 1 : size(A, 2)
                iprop = full(A(:, ind));
                iprop = unique(iprop(iprop>0));
                if ~isempty(iprop)
                    assert(numel(iprop) == 1, 'There should be only one value for a given row');
                    propfunction = cgt.model.propertyFunctionList{iprop};
                    fncallstr = propfunction.getFunctionCallString();
                    if ~isempty(fncallstr)
                        funcCallList{end + 1} = fncallstr;
                    end
                end
            end

            % we remove duplicates (the same function can be used for several output variables)
            [~, ia, ic] = unique(funcCallList, 'first');

            ia = sort(ia); % Not sure if it is necessary, but in this way we make sure the ordering for which  funcCallList was built is respected
            funcCallList = funcCallList(ia);

        end

        function printOrderedFunctionCallList(cgt)
        % Print the function calls ordered in the right order.

            funcCallList = cgt.getOrderedFunctionCallList();

            fprintf('Function call list\n');
            for ind = 1 : numel(funcCallList)
                fprintf('%s\n', funcCallList{ind});
            end

        end

        function printSubModelNames(cgt, model, parents)

            if nargin < 2
                model = cgt.model;
                parents = {};
            end

            modelnames = model.getSubModelNames();
            for imodelname = 1 : numel(modelnames)
                modelname = modelnames{imodelname};

                subparents = horzcat(parents, {modelname});
                fprintf('%s\n', strjoin(subparents, '.'))
                submodel = model.(modelname);

                cgt.printSubModelNames(submodel, subparents);

            end

        end
        
        function help(cgt, varargin)
        % print help to terminal to get an overview of all the interactive functions

            cgt = cgt.setupFuncDocs();
            
            functionDocs = cgt.functionDocs;
            names = cellfun(@(functionDoc) functionDoc.name, functionDocs, 'un', false);

            
            if nargin > 1

                option = varargin{1};
                
                parsed = false;

                if ismember(option, {'printAll', 'oneline', 'help', 'printFunctions'})
                    parsed = true;
                    if strcmp(option, 'oneline')
                        oneline = true;
                        option = 'printAll';
                    end
                end

                if ~parsed

                    funcinds = regexpSelect(names, option);
                    if  ~isempty(inds)
                        option = 'printSelected'
                        parsed = true;
                    end

                end

                assert(parsed, 'option could not be parsed');

            else
                option = 'printAll';
            end

            if nargin > 1 && strcmp(varargin{end}, 'oneline')
                oneline = true;
            else
                oneline = false;
            end
            
            parfill = ParagraphFiller('parlength', 80);
            if oneline
                parfill.parlength = Inf;
            end
            
            if ismember(option, {'printAll', 'help'})

                fprintf('Help for interactive use of the computational graph tool\n\n');

                str = 'The computational graph can be used interactively to discover and help the design of new models.';
                parfill.print(str);
                fprintf('\n');

                str = 'The help function can take the following arguments:';
                parfill.print(str);
                fprintf('\n');

                args = {};

                arg.name = 'help';
                arg.str = 'print only instruction for the help function';
                args{end + 1} = arg;
                
                arg.name = 'printAll';
                arg.str = 'print everything, particular help for all of the interactive functions. This is the default argument.';
                args{end + 1} = arg;                
                
                arg.name = 'interactive_function';
                arg.str = 'print the help for te given interactive function. A substring can be given resulting in printing all the functions that match this substring.';
                args{end + 1} = arg;

                largs    = cellfun(@(arg) strlength(arg.name), args);
                maxlargs = max(largs) + 2;
                
                for iarg = 1 : numel(largs)

                    arg = args{iarg};
                    
                    lines = parfill.getLines(arg.str);

                    formatstr = sprintf('%%-%ds %%s\n', maxlargs);
                    fprintf(formatstr, ['''', arg.name, ''''], lines{1});
                    for iline = 2 : numel(lines)
                        fprintf(formatstr, '', lines{iline});
                    end
                    if ~oneline
                        fprintf('\n');
                    end

                end

                fprintf('\n');
                str = 'if the string ''online'' is added at the end of the argument list, the output will not be formated but written as a single line (shorter output)';
                parfill.print(str);
                
            end

            if ismember(option, {'printAll', 'printFunctions'})
                
                option   = 'printSelected';
                funcinds = (1 : numel(functionDocs));
                
            end

            if strcmp(option, 'printSelected')

                fprintf('\n')
                parfill.print('Description of interactive functions for  the computational graph');
                fprintf('\n')                
                
                functionDocs = functionDocs(funcinds);

                callstrs  = cellfun(@(functionDoc) functionDoc.callstr, functionDocs, 'un', false);
                lcallstrs = cellfun(@(callstr) strlength(callstr), callstrs(funcinds));
                maxl      = max(lcallstrs);

                if ~oneline
                    parfill.parlength = 60;
                end

                for ifunc = 1 : numel(functionDocs)

                    functionDoc = functionDocs{ifunc};
                    
                    lines = parfill.getLines(functionDoc.docstring);

                    formatstr = sprintf('%%-%ds %%s\n', maxl);
                    fprintf(formatstr, functionDoc.callstr, lines{1});
                    for iline = 2 : numel(lines)
                        fprintf(formatstr, '', lines{iline});
                    end
                    if ~oneline
                        fprintf('\n');
                    end

                end
            end
            
        end

        function cgt = setupFuncDocs(cgt)

            functionDocs = {};

            % printVarNames

            docstring = 'This function lists the name of all the variables declared in the model. They corresponds to the name of the nodes in the computational graph. When the function is called with an argument, it select the variables whose name is matched by the argument, in the sense that the argument is a substring of the variable name';

            functionDoc.name      = 'printVarNames';
            functionDoc.callstr   = 'cgt.printVarNames';
            functionDoc.docstring = docstring;

            functionDocs{end + 1} = functionDoc;

            % openPropFunction(cgt, nodename)

            docstring = 'This function sends you in the matlab editor to the place in the code where the variable is updated. If the namedoes not match a unique variable name, a list of matching ones is given and the user should enter the number given in the list for the variable he/she is interested in';

            functionDoc.name      = 'openPropFunction';
            functionDoc.callstr   = 'cgt.openPropFunction(name)';
            functionDoc.docstring = docstring;

            functionDocs{end + 1} = functionDoc;

            % printRootVariables

            docstring = 'This function prints the name of the variables that are detected as roots in the graph. Those variables will correspond to the primary variables, except those that have been declared as static. The static variables are not updated by the Newton solver as they are not considered as unknown. The developper should take care of updating those explicitly';

            functionDoc.name      = 'printRootVariables';
            functionDoc.callstr   = 'cgt.printRootVariables';
            functionDoc.docstring = docstring;

            functionDocs{end + 1} = functionDoc;

            % printTailVariables

            docstring = 'This function prints the name of the variables that are detected as the tails in the graph. Those variables will correspond to the equations, except those that have been declared as extra variables. An equation variable, also called residual, is a variable that the solver will seek for its value to equal zero. We use Newton algorithm for that. The extra variables are variables that do not enter into the evaluation of the residuals but are usefull in a postprocessing of the solution.';

            functionDoc.name      = 'printTailVariables';
            functionDoc.callstr   = 'cgt.printTailVariables';
            functionDoc.docstring = docstring;

            functionDocs{end + 1} = functionDoc;

            % printChildDependencyList

            docstring = 'This functions prints the list of the function calls that will be used to update all the variables up to the residuals. The list is ordered to obey the dependency relationships that are declared in the graph';

            functionDoc.name      = 'printOrderedFunctionCallList';
            functionDoc.callstr   = 'cgt.printOrderedFunctionCallList';
            functionDoc.docstring = docstring;

            functionDocs{end + 1} = functionDoc;
            
            % printChildDependencyList

            docstring = 'Given a variable, prints the list of the variables that depends on it (children in the directed graph). The variables are given with the distance to the input variable. The distance gives an idea on how far the variable is in the evaluation tree. More precisely, in an acyclic directed graph, the distance corresponds to number of nodes that separates two nodes using the shortest path to connect them.';

            functionDoc.name      = 'printChildDependencyList';
            functionDoc.callstr   = 'cgt.printChildDependencyList(name)';
            functionDoc.docstring = docstring;

            functionDocs{end + 1} = functionDoc;
            
            % printParentDependencyList

            docstring = 'Given a variable, prints the list of the variables that the variables depends on (parents in the directed graph). The variables are given with the distance to the input variable. The distance gives an idea on how far the variable is in the evaluation tree. More precisely, in an acyclic directed graph, the distance corresponds to number of nodes that separates two nodes using the shortest path to connect them.';

            functionDoc.name      = 'printParentDependencyList';
            functionDoc.callstr   = 'cgt.printParentDependencyList(name)';
            functionDoc.docstring = docstring;

            functionDocs{end + 1} = functionDoc;

            % printPropFunctionCallList

            docstring = 'Given a variable name, this function prints the list of the function calls that will be used to update this variables. The list is ordered to obey the dependency relationships that are declared in the graph';

            functionDoc.name      = 'printPropFunctionCallList';
            functionDoc.callstr   = 'cgt.printPropFunctionCallList(name)';
            functionDoc.docstring = docstring;

            functionDocs{end + 1} = functionDoc;

            % printPropFunction

            docstring = 'Given a variable name, this function prints the correponding PropFunction object attached to this variable. It includes the function call and the list of the arguments the function depends on';

            functionDoc.name      = 'printPropFunction';
            functionDoc.callstr   = 'cgt.printPropFunction(name)';
            functionDoc.docstring = docstring;

            functionDocs{end + 1} = functionDoc;

            % printSubModelNames

            docstring = 'This function prints the list of the sub-models that constitutes the main model';

            functionDoc.name      = 'printSubModelNames';
            functionDoc.callstr   = 'cgt.printSubModelNames';
            functionDoc.docstring = docstring;

            functionDocs{end + 1} = functionDoc;

            cgt.functionDocs = functionDocs;
        end

        function exportDOT(cgt, filename)

            function str = getnodelabel(nodename)
                strs = split(nodename, '.');
                str = join(strs, '\n');
                str = str{1};
            end
            
            fid = fopen(filename, 'w');
            fprintf(fid, 'digraph G {\n');
            for ind = 1 : numel(cgt.nodenames)
                fprintf(fid, '%d [label = "%s"]\n', ind, getnodelabel(cgt.nodenames{ind}));
            end
            [ii, jj, ss] = find(cgt.adjencyMatrix);

            propfuncs = cgt.model.propertyFunctionList;
            function str = getproplabel(indprop)
                propfunc = cgt.model.propertyFunctionList{indprop};
                str = func2str(propfunc.fn);
            end
            
            for ind = 1 : numel(ii)
                fprintf(fid, '%d -> %d [label = "%s"]\n', ii(ind), jj(ind), getproplabel(ss(ind)));
            end
            
            fprintf(fid, '}\n');
            fclose(fid);
            
        end
        
        function exportCytoscape(cgt, filename)

            function str = getnodelabel(nodename)
                strs = split(nodename, '.');
                str = join(strs, ' ');
                str = ['"', str{1}, '"'];
            end

            A = cgt.adjencyMatrix;
            nodenames = cgt.nodenames;

            %% print root variables after removing the variables that were declared as static in the model

            rootnames = nodenames(all(A == 0, 1));
            staticnames = cgt.getStaticVarNames();

            ind = ismember(rootnames, staticnames);
            rootind = find(~ind);
            rootnames = rootnames(rootind);
            isroot = ismember(nodenames, rootnames);

            extravarnameinds = cgt.extraVarNameInds;

            tailinds = find(all(A' == 0, 1));
            isadded  = ismember(tailinds, extravarnameinds);
            tailinds = tailinds(~isadded);
            istail = false(numel(nodenames), 1);
            istail(tailinds) = true;
            
            nodefilename = [filename, '_nodes.csv'];
            
            nodefid = fopen(nodefilename, 'w');
            fprintf(nodefid, 'node id, node name,color, border width\n');
            for ind = 1 : numel(cgt.nodenames)
                if isroot(ind)
                    % hard-coded for the moment, just for testing ...
                    color = '#0B0A0A';
                    borderwidth = 100;
                elseif istail(ind)
                    color = '#00FF22';
                    borderwidth = 100;                    
                else
                    color = '#F7E7E7';
                    borderwidth = 2;
                end
                fprintf(nodefid, '%d,%s,%s,%d\n', ind, getnodelabel(cgt.nodenames{ind}), color, borderwidth);
            end
            fclose(nodefid);
            
            siffilename  = [filename, '.sif'];
            edgefilename = [filename, '_edgelabels.csv'];

            siffid  = fopen(siffilename, 'w');
            edgefid = fopen(edgefilename, 'w');
            
            fprintf(edgefid, 'edge id,edge name\n');

            [ii, jj, ss] = find(cgt.adjencyMatrix);

            function str = getproplabel(indprop)
                propfunc = cgt.model.propertyFunctionList{indprop};
                str = func2str(propfunc.fn);
            end
            
            for ind = 1 : numel(ii)
                
                fprintf(siffid, '%d interaction %d\n', ii(ind), jj(ind));
                fprintf(edgefid, '%d,%s\n', ind, getproplabel(ss(ind)));
                
            end

            fclose(siffid);
            fclose(edgefid);
            
        end
        
    end

    methods (Static)

        function state = setupState(state, modelspace)

            if isempty(modelspace)
                state = [];
            else
                state.(modelspace{1}) = ComputationalGraphTool.setupState(state, modelspace(2 : end));
            end
            
        end
            
        function nodenames = getNodeName(varname)
        % convert variable names to graph node names
        % TODO : we should enforce that the same function is used in setupGraph (important!)
            nodenames = {};
            varname_s = varname.resolveIndex();
            for ind = 1 : numel(varname_s)
                nodenames{end + 1} = varname_s{ind}.getIndexedFieldname();
            end

        end

        function nodenames = getNodeNames(varnames)
        % convert variable names to graph node names
        % TODO : we should enforce that the same function is used in setupGraph (important!)
            nodenames = {};
            for ind = 1 : numel(varnames)
                varname = varnames{ind};
                nodenames = horzcat(nodenames, ...
                                    ComputationalGraphTool.getNodeName(varname));
            end

        end

        function printHeader(headertxt, n)
        % Minor utility function used in this class to print a header with a number
            str = sprintf('\n%d %s', n, headertxt);
            fprintf('%s:\n', str);
            fprintf('%s\n', repmat('-', length(str), 1));
        end
    end

end




%{
Copyright 2021-2023 SINTEF Industry, Sustainable Energy Technology
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
