classdef BatteryInputParams1D < BatteryInputParams
    
    properties
        % we add those here because it facilitates setup
        sepnx  
        nenx   
        penx   
        ccnenx 
        ccpenx 
    end
    
    methods
        
        function params = BatteryInputParams1D()
            
            fac = 1;
            params.sepnx  = 30*fac;
            params.nenx   = 30*fac;
            params.penx   = 30*fac;
            params.ccnenx = 20*fac;
            params.ccpenx = 20*fac;
            
            params = params@BatteryInputParams();
        
        end
        
        function params = setupVariousParams(params)
            
            params.SOC  = 0.5;
            params.T    = 298.15;
            params.J    = 0.1;
            params.Ucut = 2;
        end

        function params = setupSubModels(params)
            
            sepnx  = params.sepnx;
            nenx   = params.nenx;
            penx   = params.penx;
            ccnenx = params.ccnenx;
            ccpenx = params.ccpenx;

            nxs = [ccnenx; nenx; sepnx; penx; ccpenx];

            xlength = 1e-6*[10; 100; 50; 80; 10];
            ylength = 1e-2;

            x = xlength./nxs;
            x = rldecode(x, nxs);
            x = [0; cumsum(x)];

            G = tensorGrid(x);
            G = computeGeometry(G); 
            params.G = G;

            %% setup elyte
            nx = sum(nxs);

            istart = ccnenx;
            ncells = nenx + sepnx + penx;
            cells = istart + (1 : ncells)';
            params.elyte = orgLiPF6('elyte', G, cells);

            %% setup ne
            istart = ccnenx;
            ncells = nenx;
            cells = istart + (1 : ncells)';
            params.ne = graphiteElectrode('ne', G, cells);

            %% setup pe
            istart = ccnenx + nenx + sepnx + 1;
            ncells = penx;
            cells = istart + (1 : ncells)';
            params.pe = nmc111Electrode('pe', G, cells);

            %% setup ccne
            istart = 1;
            ncells = ccnenx;
            cells = istart + (1 : ncells)';
            params.ccne = currentCollector('ccne', G, cells);

            %% setup ccpe
            istart = ccnenx + nenx + sepnx + penx + 1;
            ncells = ccpenx;
            cells = istart + (1 : ncells)';
            params.ccpe = currentCollector('ccpe', G, cells);

            %% setup sep
            istart = ccnenx + nenx + 1;
            ncells = sepnx;
            cells = istart + (1 : ncells)';
            params.sep = celgard2500('sep', G, cells);
        
        end
        
        function coupTerm = setupCcneBcCoupTerm(params)

            compnames = {'ccne'};
            coupTerm = couplingTerm('bc-ccne', compnames);
            coupTerm.couplingfaces = 1;
            coupTerm.couplingcells = 1;

        end

        function coupTerm = setupCcneNeCoupTerm(params)

            ccnenx = params.ccnenx;
            
            compnames = {'ccne', 'ne'};
            coupTerm = couplingTerm('ccne-ne', compnames);
            coupTerm.couplingfaces = [ccnenx, 1];
            coupTerm.couplingcells = [ccnenx, 1];

        end
        
        function coupTerm = setupCcpeBcCoupTerm(params)

            ccpenx = params.ccpenx;

            compnames = {'ccpe'};
            coupTerm = couplingTerm('bc-ccpe', compnames);
            coupTerm.couplingfaces = ccpenxx;
            coupTerm.couplingcells = ccpenxx;

        end

        function coupTerm = setupCcpePeCoupTerm(params)

            penx = params.penx;
            
            compnames = {'ccpe', 'pe'};
            coupTerm = couplingTerm('ccpe-pe', compnames);
            coupTerm.couplingfaces = [1, penx];
            coupTerm.couplingcells = [1, penx];
            
        end

        function coupTerm = setupNeElyteCoupTerm(params)
            
            nenx   = params.nenx;

            compnames = {'ne', 'elyte'};
            coupTerm = couplingTerm('ne-elyte', compnames);
            cells1 = (1 : nenx)';
            cells2 = (1 : nenx)';
            coupTerm.couplingcells =  [cells1, cells2];
            coupTerm.couplingfaces = []; % no coupling throug faces. We set it as empty
            
        end
        
        function coupTerm = setupPeElyteCoupTerm(params)
            
            sepnx  = params.sepnx;
            nenx   = params.nenx;
            penx   = params.penx;
            
            compnames = {'pe', 'elyte'};
            coupTerm = couplingTerm('pe-elyte', compnames);
            cells1 = (1 : penx)';
            cells2 = nenx + sepnx + (1 : penx)';
            coupTerm.couplingcells = [cells1, cells2];
            coupTerm.couplingfaces = []; % no coupling between faces
            
        end

    end
    
end
