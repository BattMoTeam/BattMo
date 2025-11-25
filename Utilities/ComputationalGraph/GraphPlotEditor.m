classdef GraphPlotEditor

    properties
        h % graphplot handle
    end
    
    methods

        function gpe = GraphPlotEditor(h)
            
            gpe.h = h;
            
        end

        function inds = getNodeIndex(gpe, str)

            nl = gpe.h.NodeLabel;
            if iscell(str)
                
                inds = [];
                for istr = 1 : numel(str)
                    inds = [inds, gpe.getNodeIndex(str{istr})];
                end
                inds = unique(inds);
                
            else
                inds = regexpSelect(nl, str);
            end
            
        end

        function printNodes(gpe, str)

            inds = gpe.getNodeIndex(str);
            for iind  = 1 : numel(inds)
                fprintf('%d : %s\n', inds(iind), gpe.h.NodeLabel{inds(iind)});
            end
            
        end

        function alignNodes(gpe, nodes, alignnode, direction)

            inds = gpe.getNodeIndex(nodes);
            if isnumeric(alignnode)
                alignind = gpe.getNodeIndex(nodes{alignnode});
            else
                alignind = gpe.getNodeIndex(alignnode);
            end

            assert(numel(alignind) == 1, 'align node should be unique')

            switch direction
              case 'vert'
                ref = gpe.h.XData(alignind);
              case 'horz'
                ref = gpe.h.YData(alignind);
              otherwise
                error('direction not recognized');
            end

            for iind = 1 : numel(inds)
                ind = inds(iind);
                switch direction
                  case 'vert'
                    gpe.h.XData(ind) = ref;
                  case 'horz'
                    gpe.h.YData(ind) = ref;
                  otherwise
                    error('direction not recognized');
                end
            end
            
        end
        
        function moveNode(gpe, str, dx)
            
            h = gpe.h;
            inds = gpe.getNodeIndex(str);
            
            if numel(inds) > 1
                for iind  = 1 : numel(inds)
                    fprintf('%s\n', h.NodeLabel{inds(iind)});
                end
                fprintf('specify node\n')
            else
                h.XData(ind) = h.XData(ind) + dx(1);
                h.YData(ind) = h.YData(ind) + dx(2);
            end
            
        end

        function moveNodes(gpe, inds, dx)

            h = gpe.h;

            for iind = 1 : numel(inds)
                
                ind = inds(iind);
                fprintf('%s\n', h.NodeLabel{ind});
                
                h.XData(ind) = h.XData(ind) + dx(1);
                h.YData(ind) = h.YData(ind) + dx(2);
                
            end
            
        end
        
        function moveNodesByNames(gpe, str, dx)
            
            h = gpe.h;
            inds = gpe.getNodeIndex(str);
            gpe.moveNodes(inds, dx);

        end
    end
    
end

