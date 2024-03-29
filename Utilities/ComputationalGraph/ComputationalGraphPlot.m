classdef ComputationalGraphPlot < handle

    properties (SetAccess = private)

        computationalGraphTool
        nodenames % local copy from computationalGraphTool
        A         % local copy from computationalGraphTool
        
    end
    
    properties

        filters   % list of filters

        nodeinds % current index
        
        markStaticVariables % default = true
        plotOptions
        
    end
    
    methods
        
        function cgp = ComputationalGraphPlot(computationalGraphTool)

            cgt = computationalGraphTool;
            
            cgp.computationalGraphTool = cgt;

            cgp.A         = cgt.adjencyMatrix;
            cgp.nodenames = cgt.nodenames;

            cgp.nodeinds = (1 : numel(cgp.nodenames))';
            
            cgp.markStaticVariables = true;
            cgp.plotOptions         = {};
            
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
                    str = sprintf('%s: ''%s''', filter{1}, filter{2});
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

            nodenames = cgp.nodenames;
            A         = cgp.A;


            if iscell(filter)
                action = filter{1};
                regstr = filter{2};
            else
                action = filter;
            end
            
            switch action

              case 'select'
                
                inds = regexpSelect(nodenames(cgp.nodeinds), regstr);
                cgp.nodeinds = cgp.nodeinds(inds);
                
              case 'remove'
                
                inds = regexpSelect(nodenames(cgp.nodeinds), regstr);
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
