classdef LengthSetter1D

    properties
        
        startinds % cell array with node indices of start of each component
        endinds   % cell array with node indices of end of each component
        nnodes    % total number of nodes
    end
    
    methods
        
        function lsr = LengthSetter1D(gridGenerator)
        %  gridGenerator is the BatteryGenerator1D structure
            assert(~gridGenerator.include_current_collectors, 'we do not include current collectors')

            gen = gridGenerator; % shortcut
            
            startinds{1} = 1;
            endinds{1}   = startinds{1} + gen.nenx;
            startinds{2} = endinds{1};
            endinds{2}   = startinds{2} + gen.sepnx;
            startinds{3} = endinds{2};
            endinds{3}   = startinds{3} + gen.penx;

            lsr.startinds = startinds;
            lsr.endinds   = endinds;
            lsr.nnodes    = endinds{3};
            
        end
        
        function model = setLength(lsr, model, v, comp)
        % v is length value and comp is one component name

            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            sep = 'Separator';

            ind = ismember({ne, sep, pe}, comp);

            startind = lsr.startinds{ind};
            endind   = lsr.endinds{ind};
            nnodes   = lsr.nnodes;
            
            x = model.G.parentGrid.tPFVgeometry.nodes.coords;
            x = x + 0*v.*x; % Just to make sure x gets AD if v is AD
            
            coef = v./(x(endind) - x(startind));    % multiplier coefficient for the component nodes, see below
            delta = v - (x(endind) - x(startind)); % difference to add for the following components 
            
            ind = (startind : endind);
            x(ind) = x(startind) + coef.*(x(ind) - x(startind));
            
            ind = ((endind + 1) : nnodes);
            if any(ind)
                x(ind) = x(ind) + delta;
            end
            
            model.G.parentGrid.tPFVgeometry.nodes.coords = x;
            
            model.G.parentGrid.updateTPFgeometry();
            
        end
        
    end
    
end
