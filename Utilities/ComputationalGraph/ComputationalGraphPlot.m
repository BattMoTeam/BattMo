classdef ComputationalGraphPlot < handle

    properties (SetAccess = private)

        computationalGraphTool
        nodenames % names of the nodes, copy from computationalGraphTool, with added comment for the static variable
        A         % local copy from computationalGraphTool
        
    end
    
    properties

        filters   % list of filters

        nodeinds % current index
        
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
                    varnameind = cgt.findVarName(varname);
                    [varnameinds, ~, ~, propdeplevels, ~, rootdeplevels] = getDependencyVarNameInds(varnameind, B);
                    levels = [rootdeplevels; propdeplevels];  

                    varnameinds = varnameinds(levels <= level);
                    
                    varnames = union(varnames, cgp.nodenames(varnameinds));
                    
                end

                selection = {'set', varnames};
                return
                
            end

        end

        function printSelector(cgp, selector)

            indent0 = '  ';
            
            function lines = setupLines(selector, indent)

                selectortype = selector{1};

                switch selectortype

                  case 'select'
                    
                    lines = sprintf('%s%s ''%s''', indent, selectortype, selector{2});

                    return
                  case {'and', 'or'}

                    lines{1} = sprintf('%s%s', indent, selectortype);
                    subselectors = selector{2};
                    for isel = 1 : numel(subselectors)
                        lines = horzcat(lines, setupLines(subselectors{isel}, [indent, indent0]));
                    end
                    return

                  case {'parents', 'children'}

                    level = selector{2};
                    lines{1} = sprintf('%s%s (level %d)', indent, selectortype, level);
                    lines = horzcat(lines, setupLines(selector{3}, [indent, indent0]))
                    return
                end
                

            end

            lines = setupLines(selector, '');

            for iline = numel(lines) : -1 : 1  
                fprintf('%s\n', lines{iline});
            end
            
        end
        
        function cgp = addFilters(cgp, filters, varargin)

            opt = struct('doPlot', true);
            [opt, extras] = merge_options(opt, varargin{:});

            for ifilter = 1 : numel(filters)
                filter = filters{ifilter};
                cgp = cgp.addFilter(filter, 'doPlot', false, extras{:});
            end
            
            if opt.doPlot
                
                cgp.plot();
                
            end
            
            
        end
        
        function cgp = addFilter(cgp, filter, varargin)

            opt = struct('doPlot', true, ...
                         'printFilters', true);
            opt = merge_options(opt, varargin{:});

            cgp.filters{end + 1} = filter;

            if opt.doPlot
                cgp = cgp.applyFilters('doPlot', opt.doPlot, ...
                                       'printFilters', opt.printFilters);
            end
            
        end


        function cgp = printFilters(cgp)

            filters = cgp.filters;

            fprintf('\nFilter list:\n\n')
            for ifilter = 1 : numel(filters)
                filter = filters{ifilter};
                if iscell(filter)
                    if isnumeric(filter{2})
                        str = sprintf('%s: (not printed)', filter{1});
                    else
                        str = sprintf('%s: ''%s''', filter{1}, filter{2});
                    end
                else
                    str = filter;
                end
                fprintf('%d) %s\n', ifilter, str);
            end

        end
        
        function cgp = reset(cgp, filter, varargin)
            
            opt = struct('doPlot', true);
            opt = merge_options(opt, varargin{:});            
            
            cgp.nodeinds = (1 : numel(cgp.nodenames))';

            if nargin > 1
                cgp = cgp.applyFilter(filter, varargin{:});
                return
            end
            
            if opt.doPlot
                
                cgp.plot();
                
            end

        end

        function cgp = resetFilters(cgp)

            cgp = cgp.reset();
            cgp.filters = {};
            
        end

        function cgp = applyFilters(cgp, varargin)
            
            opt = struct('doPlot', true, ...
                         'printFilters', true);
            [opt, extras] = merge_options(opt, varargin{:});

            cgp = cgp.reset();
            
            filters = cgp.filters;
            
            for ifilter = 1 : numel(filters)
                filter = filters{ifilter};
                [cgp, filter] = cgp.applyFilter(filter, 'doPlot', false, extras{:});
                filters{ifilter} = filter;
            end

            cgp.filters = filters;
            
            if opt.printFilters
                cgp.printFilters();
            end

            if opt.doPlot
                cgp.plot();
            end
            
        end

        function cgp = removeLast(cgp, varargin)
            
            opt = struct('doPlot', true, ...
                         'printFilters', true);
            opt = merge_options(opt, varargin{:});

            cgp.filters = cgp.filters(1 : end - 1);

            cgp = cgp.applyFilters('doPlot', opt.doPlot, ...
                                   'printFilters', opt.printFilters);
            
        end

        function [cgp, filter] = applyFilter(cgp, filter, varargin)
            
            opt = struct('doPlot', true);
            opt = merge_options(opt, varargin{:});            

            nodenames = cgp.computationalGraphTool.nodenames;
            A         = cgp.A;


            if iscell(filter)
                action = filter{1};
                arg = filter{2};
            else
                action = filter;
            end
            
            switch action

              case 'select'
                
                inds = regexpSelect(nodenames(cgp.nodeinds), arg);
                cgp.nodeinds = cgp.nodeinds(inds);

              case 'add nodes'
                
                nodeinds = vertcat(arg, cgp.nodeinds);
                cgp.nodeinds = unique(nodeinds);

              case 'strict-select'
                
                inds = strcmp(nodenames(cgp.nodeinds), arg);
                cgp.nodeinds = cgp.nodeinds(inds);
                
              case 'remove'
                
                inds = regexpSelect(nodenames(cgp.nodeinds), arg);
                cgp.nodeinds(inds) = [];
                
              case 'addParents'

                nodeinds = getDependencyVarNameInds(cgp.nodeinds, A);
                cgp.nodeinds = unique(nodeinds);
                
              case 'addChildren'

                nodeinds = getDependencyVarNameInds(cgp.nodeinds, A');
                cgp.nodeinds = unique(nodeinds);

              otherwise

                filter = {'select', filter};
                cgp = applyFilter(cgp, filter, varargin{:});
                
            end

            if opt.doPlot
                
                cgp.plot();
                
            end
            
        end

        function help(cgp)
        % Print help for interactive use

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

        function plotRelatives(cgp, varname, relative)
            
            cgt = cgp.computationalGraphTool;
            
            cgp.resetFilters();
            varnameind = cgt.findVarName(varname);
            if isempty(varnameind)
                fprintf('no node found matching expression\n');
                return
            end
            nodename = cgt.nodenames{varnameind};

            filteropts = {'printFilters', false, 'doplot', false};

            cgp.addFilter({'strict-select', nodename}, filteropts{:});

            
            switch relative
              case 'children'
                cgp.addFilter('addChildren', filteropts{:});
              case 'parents'
                cgp.addFilter('addParents', filteropts{:});
              case 'both'
                A = cgp.A;
                nodeinds = getDependencyVarNameInds(varnameind, A);
                nodeinds = vertcat(nodeinds, getDependencyVarNameInds(varnameind, A'));
                cgp.addFilter({'add nodes', unique(nodeinds)}, filteropts{:});
              otherwise
                error('relative option not recognized');
            end

            cgp.applyFilters();
                
        end

        function plotBoth(cgp, varname)

            cgp.plotRelatives(varname, 'both');
            
        end
        
        function plotChildren(cgp, varname)

            cgp.plotRelatives(varname, 'children');

        end

        function plotParents(cgp, varname)

            cgp.plotRelatives(varname, 'parents');
                
        end
        
        function h = plot(cgp)

            nodeinds = cgp.nodeinds;

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
