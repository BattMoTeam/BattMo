classdef ComputationalGraphPlot < handle

    properties (SetAccess = private)

        computationalGraphTool
        nodenames % names of the nodes, copy from computationalGraphTool, with added comment for the static variable
        A         % local copy from computationalGraphTool
        
    end
    
    properties


        stack = {} % stack of selectors
        
        plotOptions
        
    end
    
    methods
        
        function cgp = ComputationalGraphPlot(computationalGraphTool, varargin)

            opt = struct('markStatic', true);
            opt = merge_options(opt, varargin{:});

            cgt = computationalGraphTool;
            
            cgp.computationalGraphTool = cgt;

            cgp.A = cgt.adjencyMatrix;

            nodenames = cgt.nodenames;

            if opt.markStatic
                staticprops = cgt.staticprops;
                for istat = 1 : numel(staticprops)
                    ind = staticprops{istat}.varnameind;
                    nodenames{ind} = sprintf('%s (static)', nodenames{ind});
                end
            end
            
            cgp.nodenames = nodenames;
            
            cgp.nodeinds = (1 : numel(cgp.nodenames))';
            
            cgp.plotOptions = {};
            
        end

        function booleanOperator(cgp, op, n)

            assert(ismember(op, {'and', 'or'}), 'boolean operator not recognized');

            stack = cgp.stack;

            if nargin < 3
                n = 2;
            end
            
            if numel(stack) < n
                error(sprintf('%s operation require %d elements in this call'), op, n);
            end

            cgp.stack = { {op, stack(1 : n)}, stack{n + 1 : end}};

            cgp.printStack();
            
        end

        function and(cgp, n)
            
            if nargin < 2
                n = 2;
            end

            cgp.booleanOperator('and', n);
            
        end

        function or(cgp, n)
            
            if nargin < 2
                n = 2;
            end

            cgp.booleanOperator('or', n);
            
        end

        function addFamily(cgp, branch, level)

            assert(ismember(branch, {'parents', 'children'}), 'branch should be either parents or children');

            if nargin < 3
                level = inf;
            end

            stack = cgp.stack;

            assert(numel(stack) > 0, 'stack is empty');

            cgp.stack = {{branch, level, stack{1}}, stack{2 : end}};

            cgp.printStack();
            
        end


        function parents(cgp, level)
            
            if nargin < 2
                level = inf;
            end

            cgp.addFamily('parents', level);
            
        end

        function children(cgp, level)
            
            if nargin < 3
                level = inf;
            end

            cgp.addFamily('children', level);
            
        end
        
        function select(cgp, expr)

            cgp.stack = {{'select', expr}, cgp.stack{:}};
            cgp.printStack();
            
        end

        function reset(cgp)

            cgp.stack = {};
            
        end

        function del(cgp, n)

            if nargin < 2
                n = 1;
            end

            stack = cgp.stack;

            assert(numel(stack) >= n, sprintf('I cannot remove %d elements in the stack. Stack contains %d elements', n, numel(stack)));

            cgp.stack = stack(n + 1 : end);
            
            cgp.printStack();

        end

        function delop(cgp)

            stack = cgp.stack;
            
            assert(numel(stack) >= 1, 'stack is empty');

            selector = stack{1};

            stack = stack(2 : end);

            selectortype = selector{1};

            switch selectortype

              case {'select', 'set'}

                error('no operator at bottom of the stack');

              case {'and', 'or', 'diff'}

                cgp.stack = {selector{2}{:}, stack{:}};

              case {'parents', 'children'}

                cgp.stack = {selector{3}, stack{:}};

            end

            cgp.printStack();

        end

        function dup(cgp)

            stack = cgp.stack;

            assert(numel(stack) > 0, 'stack is empty');

            cgp.stack = {stack{1}, stack{1}, stack{2 : end}};
            
            cgp.printStack();

        end
        
        
        function swap(cgp, n)

            stack = cgp.stack;

            if nargin < 2
                n = 2;
            end
            
            assert(numel(stack) >= n, 'There should be at least %d elements in the stack to swap the %dth element', n);

            inds = (1 : numel(stack));
            inds(n) = [];
            inds = [n, inds];
            
            cgp.stack = stack(inds);

            cgp.printStack();
            
        end

        function diff(cgp)

            stack = cgp.stack;

            assert(numel(stack) >= 2, 'There should be at least 2 elements in the stack to take a diff');

            cgp.stack = {{'diff', {stack{1}, stack{2}}}, stack{3 : end}};

            cgp.printStack();
            
        end
        
        
        function selection = parseSelector(cgp, selector)

            cgt = cgp.computationalGraphTool;
            
            selectiontype = selector{1};

            assert(ischar(selectiontype), 'The first element of the selector should be a string');
            
            switch selectiontype
                
              case 'set'

                return
                
              case 'select'

                nodenames = cgp.computationalGraphTool.nodenames;
                inds = regexpSelect(nodenames, selector{2});

                selection = {'set', cgp.nodenames(inds)};
                return
                
              case {'and', 'or'}
               
                varnamesets = selector{2};

                varnameset = cgp.parseSelector(varnamesets{1});
                varnames = varnameset{2};
                for iselect = 1 : numel(varnamesets)
                    varnameset = cgp.parseSelector(varnamesets{iselect});
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
                    B = cgp.A;
                  case 'children'
                    B = cgp.A';
                end
                
                level      = selector{2};
                varnameset = cgp.parseSelector(selector{3});

                varnames   = varnameset{2};

                for ivar = 1 : numel(varnames)

                    varname = varnames{ivar};

                    % setup from method getDependencyList. We could have cleaned up the implementation there.
                    varnameind = cgt.findVarName(sprintf('%s$', varname));
                    [varnameinds, ~, ~, propdeplevels, ~, rootdeplevels] = getDependencyVarNameInds(varnameind, B);
                    levels = [rootdeplevels; propdeplevels];  

                    varnameinds = varnameinds(levels <= level);
                    
                    varnames = union(varnames, cgp.nodenames(varnameinds));

                    
                end

                selection = {'set', varnames};
                return

              case 'diff'

                args = selector{2};

                assert(numel(args) == 2, 'we expect 2 arguments for a diff');

                for iarg = 1 : numel(args)
                    varnameset = cgp.parseSelector(args{iarg});
                    varnames{iarg} = varnameset{2};
                end

                varnames = setdiff(varnames{2}, varnames{1});
                selection = {'set', varnames};
                
                return
                
            end

        end

        function printStack(cgp)

            stack = cgp.stack;
            for iselector = numel(stack) : -1  : 1
                lines = cgp.setupSelectorPrint(stack{iselector});
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

        function printStackSelection(cgp)

            assert(numel(cgp.stack) > 0, 'stack is empty');
            
            cgp.printSelection(cgp.stack{1});
            
        end
        
        function printSelection(cgp, selection)

            if ~strcmp(selection{1}, 'set')
                selection = cgp.parseSelector(selection);
            end
            
            varnames = selection{2};
            fprintf('\n');
            
            for ivar = 1 : numel(varnames)
                fprintf('%s\n', varnames{ivar});
            end

        end

        function lines = setupSelectorPrint(cgp, selector)
            
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
        
        function printSelector(cgp, selector)

            lines = cgp.setupSelectorPrint(selector);
            
            for iline = numel(lines) : -1 : 1  
                fprintf('%s\n', lines{iline});
            end
            
        end

        function help(cgp)
        % Print help for interactive use

            error('update')

            str = {};
            str{end + 1} = '';
            str{end + 1} = 'Commands available for ComputationalGraphPlot object given by cgp';
            str{end + 1} = '';
            str{end + 1} = 'cgp.plot()                       : Plot graph';
            str{end + 1} = 'cgp.addFilter(regstr)            : Select the nodes that match the regular expression regstr';
            str{end + 1} = '                                   This action is added in the filter list';
            str{end + 1} = 'cgp.removeLast()                 : Remove last filter in list';
            str{end + 1} = 'cgp.resetFilters()               : Clear filter list';
            str{end + 1} = 'cgp.addFilter(''addParents'')      : Add parents of the visible nodes';
            str{end + 1} = '                                   This action is added to the filter list';
            str{end + 1} = 'cgp.addFilter(''addChildren'')     : Add children of the visible nodes';
            str{end + 1} = '                                   This action is added to the filter list';
            str{end + 1} = 'cgp.remove({''remove'', regstr})   : Remove from the visible nodes the nodes that match the regular expression';
            str{end + 1} = '                                   This action is added to the filter list';            
            str = strjoin(str, newline);

            fprintf('%s\n', str);
            
        end

        function cgt = cgt(cgp)
        % Convenience function
            cgt = cgp.computationalgraphplot;
            
        end

        
        function h = plot(cgp)

            if isempty(cgp.stack)
                selection = {'set', cgp.nodenames};
            else
                selection = cgp.parseSelector(cgp.stack{1});
            end
            
            varnames = selection{2};

            nodeinds = ismember(cgp.nodenames, varnames);
            
            A         = cgp.A(nodeinds, nodeinds);
            nodenames = cgp.nodenames(nodeinds);
            
            g = digraph(A, nodenames);
            h = plot(g, cgp.plotOptions{:});

            if nargout < 1
                clear h
            end
                
        end
        
        function h = plotModelGraph(cgp, modelname)

            cgt = cgp.computationalGraphTool;
            
            if isempty(cgt.modelnames)
                cgt = cgt.setupModelGraph();
                cgp.computationalGraphTool = cgt;
            end

            A          = cgt.modelAdjencyMatrix;
            modelnames = cgt.modelnames;

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
    end

    methods (Static)

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
