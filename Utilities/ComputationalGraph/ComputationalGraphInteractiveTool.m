classdef ComputationalGraphInteractiveTool < handle
%
% A model is essentially characterised by its computational graph. The functions to setup the graph of a model are given in the BaseModel class.
%
% The ComputationalGraphInteractiveTool is used to store the computational graph and provides utility functions to explore it in an iteractive way.
%
% From a model, you can get the computatonial graph by writing in the terminal
%
% cgti = model.cgti
%
% Then write
%
% cgti.help

% to get an overview of the different commands that are available.
% 
    properties (SetAccess = private)

        functionDocs % cell array. Each cell contains a struct describing the function documentation, with fields
                     % - name      : function name
                     % - docstring : string describing the function

    end

    properties

        computationalGraph 

        plotnodenames % names of the nodes, copy from computationalGraphTool, with added comment for the static variable
        
        stack = {} % stack of selectors (plotting tool)
        
        plotOptions

    end

    methods

        function cgti = ComputationalGraphInteractiveTool(cg, varargin)
            
            opt = struct('markStatic', true);
            opt = merge_options(opt, varargin{:});

            cgti.computationalGraph = cg;

            nodenames = cg.nodenames;

            plotnodenames = nodenames;
            if opt.markStatic
                staticprops = cg.staticprops;
                for istat = 1 : numel(staticprops)
                    ind = staticprops{istat}.varnameind;
                    plotnodenames{ind} = sprintf('%s (static)', nodenames{ind});
                end
            end
            
            cgti.plotnodenames = plotnodenames;
            
            cgti.plotOptions = {};
            
        end


        %%%%%%%%%%%%%%%%
        %% REPL tools
        %%%%%%%%%%%%%%%%
        
        function printChildDependencyList(cgti, varname)
            
            cgti.printDependencyList(varname, 'downwards');
            
        end
        
        function printParentDependencyList(cgti, varname)
            
            cgti.printDependencyList(varname, 'upwards');

        end

        
        function printDependencyList(cgti, varname, direction)
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
            

            cg = cgti.computationalGraph;
            
            varnameind = cg.findVarName(varname);
            varname = cg.varNameList{varnameind};
            
            [varnames, varnameinds, propfuncinds, distance] = cg.getDependencyList(varname, direction);

            for ivar = 1 : numel(varnameinds)
                varnameind = varnameinds(ivar);
                fprintf('%s (%d)\n', cg.nodenames{varnameind}, distance(ivar));
            end

        end


        function printPropFunctionCallList(cgti, propfunc, varargin)
        % input propfunc is either
        % - an instance of PropFunction
        % - a valid input for findPropFunction, that is, either
        %     - a VarName instance
        %     - a cell which then uses shortcuts for VarName (see implementation below)
        %     - a string giving a regexp. It will be used to select varnames by the string name
        %   In this case, findPropFunction is first run to obtain the property function. It should return a unique element, otherwise we get a warning with a list of returned elements
        % Print the list of function call (as string) that will update the property function propfunc.

            cg = cgti.computationalGraph;
            
            opt = struct('fullSignature', false);
            opt = merge_options(opt, varargin{:});

            if ~isa(propfunc, 'PropFunction')

                varname = propfunc;

                propfunc = cg.findPropFunction(varname);

                if isempty(propfunc)
                    fprintf('No property matching regexp has been found\n');
                    return
                end

                if numel(propfunc) > 1
                    fprintf('Several property functions are matching\n\n');
                    cgti.printPropFunction(varname);
                    return
                end

            end

            nodenames = cg.nodenames;

            [propfuncs, ~, varnameinds] = cg.getPropFunctionList(propfunc);

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
                        strs = horzcat(strs, cg.getNodeName(inputvarname));
                    end
                    strs = strjoin(strs, ', ');
                    varstr = sprintf('%-*s', strl + 6, nodenames{varnameind});
                    fprintf('%s [%s]\n', varstr, strs);
                end

            end

        end


        function openPropFunction(cgti, nodename)
        % Open in editor the place where the variable that matches the regexp nodename is updated. The regexp nodename
        % should return a unique match

            cg = cgti.computationalGraph;

            [propfunc, propfuncinds] = cg.findPropFunction(nodename);

            if isempty(propfunc)
                fprintf('No property matching regexp has been found\n');
                return
            end

            if numel(propfunc) > 1
                fprintf('Several property functions are matching\n\n');
                for iprop = 1 : numel(propfunc)
                    varname = propfunc{iprop}.varname;
                    nodenames = cg.getNodeName(varname);
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
            
            state         = ComputationalGraphInteractiveTool.setupState([], mn);
            state0        = state; % may be needed in case of accumulation term
            dt            = 0;     % may be needed in case of accumulation term
            drivingForces = [];    % may be needed in case of update of driving force term
            model         = cg.model; % needed in function call passed in eval
            
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

        function printVarNames(cgti, nodename)
        % Print the variable(s) that match the regexp nodename

            cg = cgti.computationalGraph;

            nodenames = cg.nodenames;

            if nargin < 2
                nodename = '.';
            end
            inds = regexpSelect(cg.nodenames, nodename);
            nodenames = nodenames(inds);
            for ind = 1 : numel(nodenames)
                fprintf('%s\n', nodenames{ind});
            end

        end


        function printPropFunction(cgti, nodename)
        % Print property function including output variable name, function name and input variable names for the variable(s)
        % that match the regexp nodename

            cg = cgti.computationalGraph;

            propfuncs = cg.findPropFunction(nodename);

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


        function printModelNames(cgti, modelname)

            cg = cgti.computationalGraph;

            if isempty(cg.modelnames)
                cg = cg.setupModelGraph();
            end

            if nargin < 2
                modelname = '.';
            end

            modelnames = cg.modelnames;
            inds = regexpSelect(modelnames, modelname);
            modelnames = modelnames(inds);

            for imodel = 1 : numel(modelnames)
                fprintf('%s\n', modelnames{imodel});
            end

        end


        function printRootVariables(cgti, nodename)
        % Print the root variables in computational graph

            cg = cgti.computationalGraph;

            A = cg.adjencyMatrix;
            nodenames = cg.nodenames;

            %% print root variables after removing the variables that were declared as static in the model

            rootnames = nodenames(all(A == 0, 1));
            staticnames = cg.getStaticVarNames();

            ind = ismember(rootnames, staticnames);
            rootinds = find(~ind);

            if nargin > 1
                inds = regexpSelect(cg.nodenames, nodename);
                rootinds = intersect(inds, rootinds);
            end
            
            cgti.printHeader('Root Variables', numel(rootinds));

            for irind = 1 : numel(rootinds)
                fprintf('%s\n', rootnames{rootinds(irind)});
            end

            if nnz(ind)
                staticinds = find(ind);
                if nargin > 1
                    inds = regexpSelect(cg.nodenames, nodename);
                    staticinds = intersect(inds, staticinds);
                end
                cgti.printHeader('Static Variables', numel(staticinds));
                for isind = 1 : numel(staticinds)
                    fprintf('%s\n', rootnames{staticinds(isind)});
                end
            end

        end

        function printTailVariables(cgti, nodename)
        % Print the tail variables of the computational graph

            cg = cgti.computationalGraph;

            A                = cg.adjencyMatrix;
            nodenames        = cg.nodenames;
            extravarnameinds = cg.extraVarNameInds;

            tailinds = find(all(A' == 0, 1));
            isadded  = ismember(tailinds, extravarnameinds);
            tailinds = tailinds(~isadded);

            if nargin > 1
                inds = regexpSelect(cg.nodenames, nodename);
                tailinds = intersect(inds, tailinds);
            end

            cgti.printHeader('Tail Variables', numel(tailinds));
            
            for itail = 1 : numel(tailinds)
                fprintf('%s\n', nodenames{tailinds(itail)});
            end

            if nargin == 1
                cgti.printHeader('Extra Variables (do not enter in residual evaluation)', numel(extravarnameinds));
                for ievar = 1 : numel(extravarnameinds)
                    fprintf('%s\n', nodenames{extravarnameinds(ievar)});
                end
            end

        end

        function printDetachedVariables(cgti)
        % Print the "detached" variables, which are the variables that are not connected to the graph. This is specially useful
        % in debugging because such variables should be eliminated from the final graph.

            cg = cgti.computationalGraph;

            A = cg.adjencyMatrix;
            nodenames = cg.nodenames;

            fprintf('Detached variables \n');
            ind1 = all(A == 0, 1);
            ind2 = all(A' == 0, 1);
            detached_nodenames = nodenames(ind1&ind2);

            for idnodename = 1 : numel(detached_nodenames)
                fprintf('%s\n', detached_nodenames(idnodename));
            end
        end


        function printOrderedFunctionCallList(cgti)
        % Print the function calls ordered in the right order.

            cg = cgti.computationalGraph;

            funcCallList = cg.getOrderedFunctionCallList();

            fprintf('Function call list\n');
            for ind = 1 : numel(funcCallList)
                fprintf('%s\n', funcCallList{ind});
            end

        end

        function printSubModelNames(cgti, model, parents)

            cg = cgti.computationalGraph;

            if nargin < 2
                model = cg.model;
                parents = {};
            end

            modelnames = model.getSubModelNames();
            for imodelname = 1 : numel(modelnames)
                modelname = modelnames{imodelname};

                subparents = horzcat(parents, {modelname});
                fprintf('%s\n', strjoin(subparents, '.'))
                submodel = model.(modelname);

                cgti.printSubModelNames(submodel, subparents);

            end

        end
        
        function help_repl(cgti, varargin)
        % print help to terminal to get an overview of all the interactive functions

            cg = cgti.computationalGraph;
            
            cgti = cgti.setupFuncDocs();
            
            functionDocs = cgti.functionDocs;
            
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

                ComputationalGraphInteractiveTool.printFunctionDocs(functionDocs, parfill, oneline);

            end
            
        end

        function cgti = setupFuncDocs(cgti)

            cg = cgti.computationalGraph;
            
            functionDocs = {};

            % printVarNames

            docstring = 'This function lists the name of all the variables declared in the model. They corresponds to the name of the nodes in the computational graph. When the function is called with an argument, it select the variables whose name is matched by the argument, in the sense that the argument is a substring of the variable name';

            functionDoc.name      = 'printVarNames';
            functionDoc.callstr   = 'cgti.printVarNames';
            functionDoc.docstring = docstring;

            functionDocs{end + 1} = functionDoc;

            % openPropFunction(cgti, nodename)

            docstring = 'This function sends you in the matlab editor to the place in the code where the variable is updated. If the namedoes not match a unique variable name, a list of matching ones is given and the user should enter the number given in the list for the variable he/she is interested in';

            functionDoc.name      = 'openPropFunction';
            functionDoc.callstr   = 'cgti.openPropFunction(name)';
            functionDoc.docstring = docstring;

            functionDocs{end + 1} = functionDoc;

            % printRootVariables

            docstring = 'This function prints the name of the variables that are detected as roots in the graph. Those variables will correspond to the primary variables, except those that have been declared as static. The static variables are not updated by the Newton solver as they are not considered as unknown. The developper should take care of updating those explicitly';

            functionDoc.name      = 'printRootVariables';
            functionDoc.callstr   = 'cgti.printRootVariables';
            functionDoc.docstring = docstring;

            functionDocs{end + 1} = functionDoc;

            % printTailVariables

            docstring = 'This function prints the name of the variables that are detected as the tails in the graph. Those variables will correspond to the equations, except those that have been declared as extra variables. An equation variable, also called residual, is a variable that the solver will seek for its value to equal zero. We use Newton algorithm for that. The extra variables are variables that do not enter into the evaluation of the residuals but are usefull in a postprocessing of the solution.';

            functionDoc.name      = 'printTailVariables';
            functionDoc.callstr   = 'cgti.printTailVariables';
            functionDoc.docstring = docstring;

            functionDocs{end + 1} = functionDoc;

            % printChildDependencyList

            docstring = 'This functions prints the list of the function calls that will be used to update all the variables up to the residuals. The list is ordered to obey the dependency relationships that are declared in the graph';

            functionDoc.name      = 'printOrderedFunctionCallList';
            functionDoc.callstr   = 'cgti.printOrderedFunctionCallList';
            functionDoc.docstring = docstring;

            functionDocs{end + 1} = functionDoc;
            
            % printChildDependencyList

            docstring = 'Given a variable, prints the list of the variables that depends on it (children in the directed graph). The variables are given with the distance to the input variable. The distance gives an idea on how far the variable is in the evaluation tree. More precisely, in an acyclic directed graph, the distance corresponds to number of nodes that separates two nodes using the shortest path to connect them.';

            functionDoc.name      = 'printChildDependencyList';
            functionDoc.callstr   = 'cgti.printChildDependencyList(name)';
            functionDoc.docstring = docstring;

            functionDocs{end + 1} = functionDoc;
            
            % printParentDependencyList

            docstring = 'Given a variable, prints the list of the variables that the variables depends on (parents in the directed graph). The variables are given with the distance to the input variable. The distance gives an idea on how far the variable is in the evaluation tree. More precisely, in an acyclic directed graph, the distance corresponds to number of nodes that separates two nodes using the shortest path to connect them.';

            functionDoc.name      = 'printParentDependencyList';
            functionDoc.callstr   = 'cgti.printParentDependencyList(name)';
            functionDoc.docstring = docstring;

            functionDocs{end + 1} = functionDoc;

            % printPropFunctionCallList

            docstring = 'Given a variable name, this function prints the list of the function calls that will be used to update this variables. The list is ordered to obey the dependency relationships that are declared in the graph';

            functionDoc.name      = 'printPropFunctionCallList';
            functionDoc.callstr   = 'cgti.printPropFunctionCallList(name)';
            functionDoc.docstring = docstring;

            functionDocs{end + 1} = functionDoc;

            % printPropFunction

            docstring = 'Given a variable name, this function prints the correponding PropFunction object attached to this variable. It includes the function call and the list of the arguments the function depends on';

            functionDoc.name      = 'printPropFunction';
            functionDoc.callstr   = 'cgti.printPropFunction(name)';
            functionDoc.docstring = docstring;

            functionDocs{end + 1} = functionDoc;

            % printSubModelNames

            docstring = 'This function prints the list of the sub-models that constitutes the main model';

            functionDoc.name      = 'printSubModelNames';
            functionDoc.callstr   = 'cgti.printSubModelNames';
            functionDoc.docstring = docstring;

            functionDocs{end + 1} = functionDoc;

            cgti.functionDocs = functionDocs;
        end

        %%%%%%%%%%%%%%%
        %% Plot tools
        %%%%%%%%%%%%%%%

        function booleanOperator(cgti, op, n)

            assert(ismember(op, {'and', 'or'}), 'boolean operator not recognized');

            stack = cgti.stack;

            if nargin < 3
                n = 2;
            end
            
            if numel(stack) < n
                error(sprintf('%s operation require %d elements in this call'), op, n);
            end

            cgti.stack = { {op, stack(1 : n)}, stack{n + 1 : end}};

            cgti.printStack();
            
        end

        function and(cgti, n)
            
            if nargin < 2
                n = 2;
            end

            cgti.booleanOperator('and', n);
            
        end

        function or(cgti, n)
            
            if nargin < 2
                n = 2;
            end

            cgti.booleanOperator('or', n);
            
        end

        function addFamily(cgti, branch, level)

            assert(ismember(branch, {'parents', 'children'}), 'branch should be either parents or children');

            if nargin < 3
                level = inf;
            end

            stack = cgti.stack;

            assert(numel(stack) > 0, 'stack is empty');

            cgti.stack = {{branch, level, stack{1}}, stack{2 : end}};

            cgti.printStack();
            
        end


        function parents(cgti, level)
            
            if nargin < 2
                level = inf;
            end

            cgti.addFamily('parents', level);
            
        end

        function children(cgti, level)
            
            if nargin < 2
                level = inf;
            end

            cgti.addFamily('children', level);
            
        end
        
        function select(cgti, expr)

            cgti.stack = {{'select', expr}, cgti.stack{:}};
            cgti.printStack();
            
        end

        function reset(cgti)

            cgti.stack = {};
            
        end

        function del(cgti, n)

            if nargin < 2
                n = 1;
            end

            stack = cgti.stack;

            assert(numel(stack) >= n, sprintf('I cannot remove %d elements in the stack. Stack contains %d elements', n, numel(stack)));

            cgti.stack = stack(n + 1 : end);
            
            cgti.printStack();

        end

        function delop(cgti)

            stack = cgti.stack;
            
            assert(numel(stack) >= 1, 'stack is empty');

            selector = stack{1};

            stack = stack(2 : end);

            selectortype = selector{1};

            switch selectortype

              case {'select', 'set'}

                error('no operator at bottom of the stack');

              case {'and', 'or', 'diff'}

                cgti.stack = {selector{2}{:}, stack{:}};

              case {'parents', 'children'}

                cgti.stack = {selector{3}, stack{:}};

            end

            cgti.printStack();

        end

        function dup(cgti)

            stack = cgti.stack;

            assert(numel(stack) > 0, 'stack is empty');

            cgti.stack = {stack{1}, stack{1}, stack{2 : end}};
            
            cgti.printStack();

        end
        
        
        function swap(cgti, n)

            stack = cgti.stack;

            if nargin < 2
                n = 2;
            end
            
            assert(numel(stack) >= n, 'There should be at least %d elements in the stack to swap the %dth element', n);

            inds = (1 : numel(stack));
            inds(n) = [];
            inds = [n, inds];
            
            cgti.stack = stack(inds);

            cgti.printStack();
            
        end

        function diff(cgti)

            stack = cgti.stack;

            assert(numel(stack) >= 2, 'There should be at least 2 elements in the stack to take a diff');

            cgti.stack = {{'diff', {stack{1}, stack{2}}}, stack{3 : end}};

            cgti.printStack();
            
        end
        
        
        function selection = parseSelector(cgti, selector)

            cg = cgti.computationalGraph;
            
            selectiontype = selector{1};

            assert(ischar(selectiontype), 'The first element of the selector should be a string');
            
            switch selectiontype
                
              case 'set'

                return
                
              case 'select'

                nodenames = cg.nodenames;
                inds = regexpSelect(nodenames, selector{2});

                selection = {'set', cg.nodenames(inds)};
                return
                
              case {'and', 'or'}
                
                varnamesets = selector{2};

                varnameset = cgti.parseSelector(varnamesets{1});
                varnames = varnameset{2};
                for iselect = 1 : numel(varnamesets)
                    varnameset = cgti.parseSelector(varnamesets{iselect});
                    switch selectiontype
                      case 'and'
                        varnames = intersect(varnames, varnameset{2});
                      case 'or'
                        varnames = union(varnames, varnameset{2});
                    end
                end

                selection = {'set', varnames};

                return

              case {'parents', 'children'}

                switch selectiontype
                  case 'parents'
                    B = cg.adjencyMatrix;
                  case 'children'
                    B = cg.adjencyMatrix';
                end
                
                level      = selector{2};
                varnameset = cgti.parseSelector(selector{3});

                varnames   = varnameset{2};

                for ivar = 1 : numel(varnames)

                    varname = varnames{ivar};

                    % setup from method getDependencyList. We could have cleaned up the implementation there.
                    varnameind = cg.findVarName(sprintf('%s$', varname));
                    [varnameinds, ~, ~, propdeplevels, ~, rootdeplevels] = getDependencyVarNameInds(varnameind, B);
                    levels = [rootdeplevels; propdeplevels];  

                    varnameinds = varnameinds(levels <= level);
                    
                    varnames = union(varnames, cg.nodenames(varnameinds));

                    
                end

                selection = {'set', varnames};
                return

              case 'diff'

                args = selector{2};

                assert(numel(args) == 2, 'we expect 2 arguments for a diff');

                for iarg = 1 : numel(args)
                    varnameset = cgti.parseSelector(args{iarg});
                    varnames{iarg} = varnameset{2};
                end

                varnames = setdiff(varnames{2}, varnames{1});
                selection = {'set', varnames};
                
                return
                
            end

        end

        function printStack(cgti)

            stack = cgti.stack;
            for iselector = numel(stack) : -1  : 1
                lines = cgti.setupSelectorPrint(stack{iselector});
                nlines = numel(lines);
                for iline = nlines : -1 : 1
                    if iline == 1
                        start = sprintf('%2d: ', iselector);
                    else
                        start = '    ';
                    end
                    fprintf('%s%s\n', start, lines{iline});
                end
            end
            
        end

        function printStackSelection(cgti)

            assert(numel(cgti.stack) > 0, 'stack is empty');
            
            cgti.printSelection(cgti.stack{1});
            
        end
        
        function printSelection(cgti, selection)

            if ~strcmp(selection{1}, 'set')
                selection = cgti.parseSelector(selection);
            end
            
            varnames = selection{2};
            fprintf('\n');
            
            for ivar = 1 : numel(varnames)
                fprintf('%s\n', varnames{ivar});
            end

        end

        function lines = setupSelectorPrint(cgti, selector)
            
            indent0 = '  ';
            
            function lines = setupLines(selector, indent)

                selectortype = selector{1};

                switch selectortype

                  case 'select'
                    
                    lines{1} = sprintf('%s%s ''%s''', indent, selectortype, selector{2});
                    return

                  case {'and', 'or', 'diff'}

                    lines{1} = sprintf('%s%s', indent, selectortype);
                    subselectors = selector{2};
                    for isel = 1 : numel(subselectors)
                        lines = horzcat(lines, setupLines(subselectors{isel}, [indent, indent0]));
                    end
                    return

                  case {'parents', 'children'}

                    level = selector{2};
                    lines{1} = sprintf('%s%s (level %d)', indent, selectortype, level);
                    lines = horzcat(lines, setupLines(selector{3}, [indent, indent0]));
                    return

                end
                

            end

            lines = setupLines(selector, '');

        end
        
        function printSelector(cgti, selector)

            lines = cgti.setupSelectorPrint(selector);
            
            for iline = numel(lines) : -1 : 1  
                fprintf('%s\n', lines{iline});
            end
            
        end
        
        function h = plot(cgti)

            cg = cgti.computationalGraph;
            
            if isempty(cgti.stack)
                selection = {'set', cg.nodenames};
            else
                selection = cgti.parseSelector(cgti.stack{1});
            end
            
            varnames = selection{2};

            nodeinds = ismember(cg.nodenames, varnames);
            
            A         = cg.adjencyMatrix(nodeinds, nodeinds);
            nodenames = cg.nodenames(nodeinds);
            
            g = digraph(A, nodenames);
            h = plot(g, cgti.plotOptions{:});

            if nargout < 1
                clear h
            end
            
        end
        
        function h = plotModelGraph(cgti, modelname)

            cg = cgti.computationalGraph;
            cg = cg.setupModelGraph();
            
            A          = cg.modelAdjencyMatrix;
            modelnames = cg.modelnames;

            if nargin > 1
                inds = regexpSelect(modelnames, modelname);
                A = A(inds, inds);
                modelnames = modelnames(inds);
            end

            g = digraph(A, modelnames);
            
            h = plot(g);

            if nargout < 1
                clear h
            end
            
        end
        
        %%%%%%%%%%%%%%%
        %% Export tools
        %%%%%%%%%%%%%%%

        function exportDOT(cgti, filename)

            cg = cgti.computationalGraph;

            function str = getnodelabel(nodename)
                strs = split(nodename, '.');
                str = join(strs, '\n');
                str = str{1};
            end
            
            fid = fopen(filename, 'w');
            fprintf(fid, 'digraph G {\n');
            for ind = 1 : numel(cg.nodenames)
                fprintf(fid, '%d [label = "%s"]\n', ind, getnodelabel(cg.nodenames{ind}));
            end
            [ii, jj, ss] = find(cg.adjencyMatrix);

            propfuncs = cg.model.propertyFunctionList;
            function str = getproplabel(indprop)
                propfunc = cg.model.propertyFunctionList{indprop};
                str = func2str(propfunc.fn);
            end
            
            for ind = 1 : numel(ii)
                fprintf(fid, '%d -> %d [label = "%s"]\n', ii(ind), jj(ind), getproplabel(ss(ind)));
            end
            
            fprintf(fid, '}\n');
            fclose(fid);
            
        end
        
        function exportCytoscape(cgti, filename)

            cg = cgti.computationalGraph;
            
            function str = getnodelabel(nodename)
                strs = split(nodename, '.');
                str = join(strs, ' ');
                str = ['"', str{1}, '"'];
            end

            A = cg.adjencyMatrix;
            nodenames = cg.nodenames;

            %% print root variables after removing the variables that were declared as static in the model

            rootnames = nodenames(all(A == 0, 1));
            staticnames = cg.getStaticVarNames();

            ind = ismember(rootnames, staticnames);
            rootind = find(~ind);
            rootnames = rootnames(rootind);
            isroot = ismember(nodenames, rootnames);

            extravarnameinds = cg.extraVarNameInds;

            tailinds = find(all(A' == 0, 1));
            isadded  = ismember(tailinds, extravarnameinds);
            tailinds = tailinds(~isadded);
            istail = false(numel(nodenames), 1);
            istail(tailinds) = true;
            
            nodefilename = [filename, '_nodes.csv'];
            
            nodefid = fopen(nodefilename, 'w');
            fprintf(nodefid, 'node id, node name,color, border width\n');
            for ind = 1 : numel(cg.nodenames)
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
                fprintf(nodefid, '%d,%s,%s,%d\n', ind, getnodelabel(cg.nodenames{ind}), color, borderwidth);
            end
            fclose(nodefid);
            
            siffilename  = [filename, '.sif'];
            edgefilename = [filename, '_edgelabels.csv'];

            siffid  = fopen(siffilename, 'w');
            edgefid = fopen(edgefilename, 'w');
            
            fprintf(edgefid, 'edge id,edge name\n');

            [ii, jj, ss] = find(cg.adjencyMatrix);

            function str = getproplabel(indprop)
                propfunc = cg.model.propertyFunctionList{indprop};
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
                state.(modelspace{1}) = ComputationalGraphInteractiveTool.setupState(state, modelspace(2 : end));
            end
            
        end

        function printFunctionDocs(functionDocs, parfill, oneline)

            callstrs  = cellfun(@(functionDoc) functionDoc.callstr, functionDocs, 'un', false);
            lcallstrs = cellfun(@(callstr) strlength(callstr), callstrs);
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
        
        function printHeader(headertxt, n)
        % Minor utility function used in this class to print a header with a number
            str = sprintf('\n%d %s', n, headertxt);
            fprintf('%s:\n', str);
            fprintf('%s\n', repmat('-', length(str), 1));
        end

        function selection = getIntersection(selection1, selection2)
        % Get the intersection of two selections
            if iscell(selection1) && iscell(selection2)
                selection = intersect(selection1, selection2);
            elseif ischar(selection1) && ischar(selection2)
                if strcmp(selection1, selection2)
                    selection = {selection1};
                else
                    selection = {};
                end
            else
                error('Selection types do not match');
            end
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
