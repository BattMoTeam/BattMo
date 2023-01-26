classdef ComputationalGraphTool

    properties
        graph
        model
        A
        varNameList
        nodenames
        staticprops
    end
    
    methods
        
        function cgt = ComputationalGraphTool(model)
            if isempty(model.propertyFunctionList)
                model = model.registerVarAndPropfuncNames();
            end
            cgt.model = model;
            [g, staticprops, varNameList] = setupGraph(model);
            cgt.graph       = g;
            cgt.staticprops = staticprops;
            cgt.varNameList = varNameList;
            % In adjacency matrix A, 
            % - column index      : output variable index (as in model.varNameList)
            % - row index         : input variable (as in model.varNameList)
            % - coefficient value : property function index as in model.propertyFunctionList
            cgt.A         = adjacency(g, 'weighted');
            cgt.nodenames = g.Nodes.Variables;

            
        end

        function nodestaticnames = getStaticVarNames(cgt)

            varnames = cgt.model.staticVarNameList;
            nodestaticnames = cgt.getNodeNames(varnames);

            % add the variables that are updated with no input arguments
            propfunctions = cgt.model.propertyFunctionList;
            varnames = {};
            for iprop = 1 : numel(propfunctions)
                propfunction = propfunctions{iprop};
                if isempty(propfunction.inputvarnames)
                    varnames{end + 1} = propfunction.varname;
                end
            end
            
            nodestaticnames = horzcat(nodestaticnames, cgt.getNodeNames(varnames));
            
        end
        
        function propfuncs = findPropFunction(cgt, nodename)
            
            nodenames = cgt.nodenames;
            A         = cgt.A;
            model     = cgt.model;
            
            indSelectedNodenames = regexp(nodenames, nodename, 'once');
            indSelectedNodenames = cellfun(@(x) ~isempty(x), indSelectedNodenames);
            propfuncinds = A(:, indSelectedNodenames);
            propfuncinds = unique(propfuncinds(:));
            propfuncinds = propfuncinds(propfuncinds > 0); % remove the zero elements

            staticprops = cgt.staticprops;
            staticnodenames = cellfun(@(staticprop) staticprop.nodename, staticprops, 'uniformoutput', false);
            indSelectedNodenames = regexp(staticnodenames, nodename, 'once');
            indSelectedNodenames = cellfun(@(x) ~isempty(x), indSelectedNodenames);
            staticprops = staticprops(indSelectedNodenames);
            staticIndPropfunctions = cellfun(@(staticprop) staticprop.propind, staticprops);

            propfuncinds = [propfuncinds; staticIndPropfunctions'];
            
            propfuncs = model.propertyFunctionList(propfuncinds);

            if numel(propfuncs) == 1
                propfuncs = propfuncs{1};
            end
            
        end

        function propfuncs = getPropFunctionList(cgt, propfunc)

            A = cgt.A;

            varname = propfunc.varname;
            % for cell-valued variable pick-up a valid index
            if (varname.dim > 1) && (ischar(varname.index))
                varname.index = 1;
            elseif isnumeric(varname.index) && (numel(varname.index) > 1)
                varname.index = varname.index(1);
            end
            
            varnameind = cgt.getVarNameIndex(varname);

            nodeindlist = getDependencyVarNameInds(varnameind, cgt.A);

            propfuncinds = [];
            
            for inode = 1 : numel(nodeindlist)

                nodeind = nodeindlist(inode);
                propfuncind = unique(A(:, nodeind));
                propfuncind = propfuncind(propfuncind > 0);
                if ~isempty(propfuncind)
                    propfuncinds(end + 1) = propfuncind;
                end
            end

            propfuncinds = propfuncinds(end : -1 : 1); 
            
            propfuncs = cgt.model.propertyFunctionList(propfuncinds);
            
        end

        function funcCallList = setPropFunctionCallList(cgt, nodename)

            foundpropfuncs = cgt.findPropFunction(nodename);

            funcCallList = {};

            if isempty(foundpropfuncs)
                fprintf('No property matching regexp has been found\n');
                return
            end            

            if isa(foundpropfuncs, 'PropFunction')
                foundpropfuncs = {foundpropfuncs};
            end
            
            
            for ifound = 1 : numel(foundpropfuncs)

                foundpropfunc = foundpropfuncs{ifound};
                propfuncs = cgt.getPropFunctionList(foundpropfunc);
                
                propstrs = {};
                
                for iprop = 1 : numel(propfuncs)

                    propfunc = propfuncs{iprop};
                    str = propfunc.functionCallSetupFn(propfunc);
                    propstrs{end + 1} =  str;
                    
                end

                [~, ia, ic] = unique(propstrs, 'first');
                ia = sort(ia);
                propstrs = propstrs(ia);

                funcCallList = horzcat(funcCallList, propstrs);
            end

        end

        function printPropFunctionCallList(cgt, nodename)

            funcCallList = cgt.setPropFunctionCallList(nodename);
            for ifunc = 1 : numel(funcCallList)
                funcCall = funcCallList{ifunc};
                fprintf('%s\n', funcCall);
            end
            
        end

        
        
        function openPropFunction(cgt, nodename)
            propfunc = cgt.findPropFunction(nodename);

            if isempty(propfunc)
                fprintf('No property matching regexp has been found\n');
                return
            end
            
            if numel(propfunc) > 1
                fprintf('Several property functions are matching\n\n');
                cgt.printPropFunction(nodename);
                return
            end

            fn = propfunc.fn;
            mn = propfunc.modelnamespace;
            mn = join(mn, '.');
            fnname = func2str(fn);
            fnname = regexp(fnname, "\.(.*)", 'tokens');
            fnname = fnname{1}{1};
            fnname = horzcat(mn, {fnname});
            fnname = join(fnname, '.');
            fnname = fnname{1};
            fnname = sprintf('cgt.model.%s', fnname);

            try
                eval(fnname)
            catch ME
                fprintf('%s\n', fnname);
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

            if nargin < 2
                nodename = '.'; % will match every thing
            end
            nodenames = cgt.nodenames;
            indSelectedNodenames = regexp(nodenames, nodename, 'once');
            indSelectedNodenames = cellfun(@(x) ~isempty(x), indSelectedNodenames);
            nodenames = nodenames(indSelectedNodenames);
            for ind = 1 : numel(nodenames)
                fprintf('%s\n', nodenames{ind});
            end
            
        end
        
        
        function printPropFunction(cgt, nodename)
            propfuncs = cgt.findPropFunction(nodename);

            if numel(propfuncs) == 1
                propfuncs = {propfuncs};
            end

            for iprop = 1 : numel(propfuncs)
                propfunc = propfuncs{iprop};
                varname = propfunc.varname;
                varname_s = varname.resolveIndex();
                outputvarstrs = {};
                for ind = 1 : numel(varname_s)
                    fullname = varname_s{ind}.getIndexedFieldname();
                    outputvarstrs{end + 1} = sprintf('state.%s', fullname);
                end
                outputvarstr = join(outputvarstrs, {', '});
                fprintf('%s <-', outputvarstr{1});
                fprintf(' (%s) ', propfuncs{iprop}.getFunctionCallString());
                inputvarnames = propfunc.inputvarnames;
                if isempty(inputvarnames)
                    fprintf('[no state field is used]\n');
                else
                    fprintf(' <- ', propfuncs{iprop}.getFunctionCallString());
                    inputvarstrs = {};
                    for ind = 1 : numel(inputvarnames)
                        varname = inputvarnames{ind};
                        inputvarstrs{end + 1} = sprintf('state.%s', varname.getFieldname());
                    end
                    inputvarstr = join(inputvarstrs, {', '});
                    fprintf('[%s]\n', inputvarstr{1});
                end
            end
            
        end

        function varnameind = getVarNameIndex(cgt, varname)

            for varnameind = 1 : numel(cgt.varNameList)
                if varname.compareVarName(cgt.varNameList{varnameind})
                    return
                end
            end

            fprintf('Variable not found\n');
            varnameind = [];

            
        end
            
        function [g, edgelabels] = getComputationalGraph(cgt, varargin)
            
            opt = struct('type'            , 'ascendant', ...
                         'oneParentOnly'   , false      , ...
                         'markStaticVars'  , true       , ...
                         'doplot'          , false      , ...
                         'includeNodeNames', []         , ...
                         'excludeNodeNames', []);
            opt = merge_options(opt, varargin{:});

            propfunctionlist = cgt.model.propertyFunctionList;
            nodenames        = cgt.nodenames;
            includeNodeNames = opt.includeNodeNames;
            excludeNodeNames = opt.excludeNodeNames;

            A = (cgt.A);
            if strcmp(opt.type, 'descendant')
                A = A';
            end
            
            if ~isempty(includeNodeNames)
                
                nodes = getDependencyVarNameIndsByName(includeNodeNames, nodenames, A, 'oneParentOnly', opt.oneParentOnly);

                nodenames = nodenames(nodes);
                A = A(nodes, nodes);
                
                if ~isempty(excludeNodeNames)
                    removenodes = regexp(nodenames, excludeNodeNames, 'once');
                    removenodes = cellfun(@(x) ~isempty(x), removenodes);
                    nodenames = nodenames(~removenodes);
                    A = A(~removenodes, ~removenodes);
                end

            end
            
            if strcmp(opt.type, 'descendant')
                A = A';
            end
            
            if opt.markStaticVars
                
                staticnames = cgt.getStaticVarNames();
                ind = ismember(nodenames, staticnames);

                if nnz(ind)
                    nodenames(ind) = cellfun(@(str) sprintf('%s (static)', str), nodenames(ind), 'un', false);
                end
                
            end
            
            g = digraph(A, nodenames);


            propinds = g.Edges.Weight;

            if nargout > 1
                edgelabels = {};
                for iprop = 1 : numel(propinds)
                    propind = propinds(iprop);
                    edgelabels{end + 1} = func2str(propfunctionlist{propind}.fn);
                end
            end

            if opt.doplot
                figure
                h = plot(g, 'nodefontsize', 10);
            end
            
        end
        
        function printRootVariables(cgt)

            A = cgt.A;
            nodenames = cgt.nodenames;
            
            %% print root variables after removing the variables that were declared as static in the model

            rootnames = nodenames(all(A == 0, 1));
            staticnames = cgt.getStaticVarNames();

            fprintf('Root variables \n');

            ind = ismember(rootnames, staticnames);
            rootnames(~ind)

            if nnz(ind)
                fprintf('static variables \n');
                rootnames(ind)
            end

        end


        function primvarnames = getPrimaryVariables(cgt)

            A = cgt.A;
            nodenames = cgt.nodenames;

            ind = find(all(A == 0, 1));
            
            rootnames = nodenames(ind);
            staticnames = cgt.getStaticVarNames();
            indstatic = ismember(rootnames, staticnames);
            ind = ind(~indstatic);
            primVarNameList = cgt.varNameList(ind);
            
            primvarnames = cellfun(@(varname) varname.getPropName(), primVarNameList, 'uniformoutput', false);

        end
        
        function printTailVariables(cgt)
            
            A = cgt.A;
            nodenames = cgt.nodenames;

            %% print tail variables
            fprintf('Tail variables \n');
            nodenames(all(A' == 0, 1))
            
        end        


        function printDetachedVariables(cgt)
        % for debugging : find out variables that are declared but not attached to any function (as argument or input)
            A = cgt.A;
            nodenames = cgt.nodenames;
            
            fprintf('Detached variables \n');
            ind1 = all(A == 0, 1);
            ind2 = all(A' == 0, 1);
            nodenames(ind1&ind2)
            
        end

        
        function funcCallList = setOrderedFunctionCallList(cgt)

            A = cgt.A;
            staticprops = cgt.staticprops;
            
            try
                p = topological_order(A);
            catch
                fprintf('You need to install matlab BGL\n');
                return
            end

            if isempty(p)
                fprintf('The graph contains cycles. It implies that some variables cannot be evaluated.\n')
            end
            
            funcCallList = {};

            if ~isempty(staticprops)
                
                for istatic = 1 : numel(staticprops)

                    iprop = staticprops{istatic}.propind;
                    propfunction = cgt.model.propertyFunctionList{iprop};
                    funcCallList{end + 1} = propfunction.getFunctionCallString();
                    
                end
                
            end
            
            for ind = 1 : numel(p)
                iprop = full(A(:, p(ind)));
                iprop = unique(iprop(iprop>0));
                if ~isempty(iprop)
                    assert(numel(iprop) == 1, 'There should be only one value for a given row');
                    propfunction = cgt.model.propertyFunctionList{iprop};
                    funcCallList{end + 1} = propfunction.getFunctionCallString();
                end
            end

            % we remove duplicates (the same function can be used for several output variables)
            [~, ia, ic] = unique(funcCallList, 'first');

            ia = sort(ia); % Not sure if it is necessary, but in this way we make sure the ordering for which  funcCallList was built is respected
            funcCallList = funcCallList(ia);


        end

        function printOrderedFunctionCallList(cgt)

            funcCallList = cgt.setOrderedFunctionCallList();
            
            fprintf('Function call list\n');
            for ind = 1 : numel(funcCallList)
                fprintf('%s\n', funcCallList{ind});
            end

        end
        
            
    end

    methods (Static)

        function nodenames = getNodeNames(varnames)
        % convert variable names to graph node names
        % TODO : we should enforce that the same function is used in setupGraph (important!)
            nodenames = {};
            for ind = 1 : numel(varnames)
                varname = varnames{ind};
                varname_s = varname.resolveIndex();
                for ind = 1 : numel(varname_s)
                    nodenames{end + 1} = varname_s{ind}.getIndexedFieldname();
                end
            end
        end
        
    end
    
    
end
    
