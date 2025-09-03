classdef Selector < handle

    properties
        
        stack = {} % stack of selectors (plotting tool)

        interactiveOptions
        
    end
    
    methods

        function slt = Selector(varargin)
            
            opt = struct('interactiveOptions', []);
            opt = merge_options(opt, varargin{:});
            
            interactiveOptions = setDefaultJsonStructField(opt.interactiveOptions, 'printStackAfterUpdate', true);

            slt.interactiveOptions = interactiveOptions;
            
        end
        
        function booleanOperator(slt, op, n)

            assert(ismember(op, {'and', 'or'}), 'boolean operator not recognized');

            stack = slt.stack;

            if nargin < 3
                n = 2;
            end
            
            if numel(stack) < n
                error(sprintf('%s operation require %d elements in this call'), op, n);
            end

            slt.stack = { {op, stack(1 : n)}, stack{n + 1 : end}};

            if slt.interactiveOptions.printStackAfterUpdate
                slt.printStack();
            end
            
        end

        function and(slt, n)
            
            if nargin < 2
                n = 2;
            end

            slt.booleanOperator('and', n);
            
        end

        function or(slt, n)
            
            if nargin < 2
                n = 2;
            end

            slt.booleanOperator('or', n);
            
        end

        function select(slt, expr)

            slt.stack = {{'select', expr}, slt.stack{:}};
            if slt.interactiveOptions.printStackAfterUpdate
                slt.printStack();
            end
            
        end

        function reset(slt)

            slt.stack = {};
            
        end

        function del(slt, n)

            if nargin < 2
                n = 1;
            end

            stack = slt.stack;

            assert(numel(stack) >= n, sprintf('I cannot remove %d elements in the stack. Stack contains %d elements', n, numel(stack)));

            slt.stack = stack(n + 1 : end);
            
            if slt.interactiveOptions.printStackAfterUpdate
                slt.printStack();
            end

        end

        function delop(slt)

            stack = slt.stack;
            
            assert(numel(stack) >= 1, 'stack is empty');

            selector = stack{1};

            stack = stack(2 : end);

            selectortype = selector{1};

            switch selectortype

              case {'select', 'set'}

                error('no operator at bottom of the stack');

              case {'and', 'or', 'diff'}

                slt.stack = {selector{2}{:}, stack{:}};

                slt.stack = {selector{3}, stack{:}};

            end

            if slt.interactiveOptions.printStackAfterUpdate
                slt.printStack();
            end

        end

        function dup(slt)

            stack = slt.stack;

            assert(numel(stack) > 0, 'stack is empty');

            slt.stack = {stack{1}, stack{1}, stack{2 : end}};
            
            if slt.interactiveOptions.printStackAfterUpdate
                slt.printStack();
            end

        end
        
        
        function swap(slt, n)

            stack = slt.stack;

            if nargin < 2
                n = 2;
            end
            
            assert(numel(stack) >= n, 'There should be at least %d elements in the stack to swap the %dth element', n);

            inds = (1 : numel(stack));
            inds(n) = [];
            inds = [n, inds];
            
            slt.stack = stack(inds);

            if slt.interactiveOptions.printStackAfterUpdate
                slt.printStack();
            end
            
        end

        function diff(slt)

            stack = slt.stack;

            assert(numel(stack) >= 2, 'There should be at least 2 elements in the stack to take a diff');

            slt.stack = {{'diff', {stack{1}, stack{2}}}, stack{3 : end}};

            if slt.interactiveOptions.printStackAfterUpdate
                slt.printStack();
            end
            
            
        end
        

        function found = find(slt, givenset, criteria)
        % Returns the indices that are matched in the given set by the criteria
            
            error('base function')
            
        end

        function str = elementToString(slt, element)
        % Returns the printed form of the element
            
            error('base function')
            
        end
        
        function selection = parseSelector(slt, givenset, selector)
            
            selectiontype = selector{1};

            assert(ischar(selectiontype), 'The first element of the selector should be a string');
            
            switch selectiontype
                
              case 'set'

                return
                
              case 'select'

                inds = slt.find(givenset, selector{2});

                selection = {'set', inds};

                return
                
              case {'and', 'or'}
                
                % list of the arguments for the boolean operator, which consist of a list of selector
                boolean_arg_selectors = selector{2};

                % We parse the first selector
                boolean_arg_set = slt.parseSelector(givenset, boolean_arg_selectors{1});

                % We recover the indices
                boolean_arg_ind = boolean_arg_set{2};
                
                for iselect = 1 : numel(boolean_arg_selectors)
                    boolean_arg_set = slt.parseSelector(givenset, boolean_arg_selectors{iselect});
                    switch selectiontype
                      case 'and'
                        boolean_arg_ind = intersect(boolean_arg_ind, boolean_arg_set{2});
                      case 'or'
                        boolean_arg_ind = union(boolean_arg_ind, boolean_arg_set{2});
                    end
                end

                selection = {'set', boolean_arg_ind};

                return

              case 'diff'

                error('not updated yet');
                
                args = selector{2};

                assert(numel(args) == 2, 'we expect 2 arguments for a diff');

                for iarg = 1 : numel(args)
                    varnameset = slt.parseSelector(args{iarg});
                    varnames{iarg} = varnameset{2};
                end

                varnames = setdiff(varnames{2}, varnames{1});
                selection = {'set', varnames};
                
                return
                
            end

        end

        function printStack(slt)

            stack = slt.stack;
            for iselector = numel(stack) : -1  : 1
                lines = slt.setupSelectorPrint(stack{iselector});
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

        function inds = parseStack(slt, givenset)

            inds = slt.parseSelector(givenset, slt.stack{1});
            
        end
        
        function subset = extract(slt, givenset)

            inds = slt.parseStack(givenset);
            subset = givenset(inds{2});
            
        end
        

        function printStackSelection(slt, givenset)

            assert(numel(slt.stack) > 0, 'stack is empty');
            
            slt.printSelection(givenset, slt.stack{1});
            
        end

        function printsel(slt, givenset)

            % shortcut
            slt.printStackSelection(givenset);
            
        end

        function printall(slt, givenset)

            selection = {'set', (1 : numel(givenset))};
            slt.printSelection(givenset, selection);
            
        end
        
        function printSelection(slt, givenset, selection)

            if ~strcmp(selection{1}, 'set')
                selection = slt.parseSelector(givenset, selection);
            end
            
            inds = selection{2};
            fprintf('\n');
            
            for ivar = 1 : numel(inds)
                ind = inds(ivar)
                fprintf('%s\n', slt.elementToString(givenset{ind}));
            end

        end

        function str = selectSelectorToString(slt, selectSelector)
        % Returns the printed form of a 'select' selector
            
            error('base function');
            
        end
        function lines = setupSelectorPrint(slt, selector)
            
            indent0 = '  ';
            
            function lines = setupLines(selector, indent)

                selectortype = selector{1};

                switch selectortype

                  case 'select'
                    
                    lines{1} = sprintf('%s%s ''%s''', indent, selectortype, slt.selectSelectorToString(selector));
                    return

                  case {'and', 'or', 'diff'}

                    lines{1} = sprintf('%s%s', indent, selectortype);
                    subselectors = selector{2};
                    for isel = 1 : numel(subselectors)
                        lines = horzcat(lines, setupLines(subselectors{isel}, [indent, indent0]));
                    end
                    return

                end
                

            end

            lines = setupLines(selector, '');

        end
        
        function printSelector(slt, selector)

            lines = slt.setupSelectorPrint(selector);
            
            for iline = numel(lines) : -1 : 1  
                fprintf('%s\n', lines{iline});
            end
            
        end


        
    end
    
end
    

