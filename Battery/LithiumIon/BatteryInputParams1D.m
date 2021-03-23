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
        
        function params = setupVariousParams(params)
            
            params.SOC  = 0.5;
            params.T    = 298.15;
            params.J    = 0.1;
            params.Ucut = 2;
            
            fac = 1;
            params.sepnx  = 30*fac;
            params.nenx   = 30*fac;
            params.penx   = 30*fac;
            params.ccnenx = 20*fac;
            params.ccpenx = 20*fac;
            
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
            params.Electrolyte = orgLiPF6('elyte', G, cells);

            %% setup ne
            istart = ccnenx;
            ncells = nenx;
            cells = istart + (1 : ncells)';
            params.NegativeElectrode = GraphiteElectrode('ne', G, cells);

            %% setup pe
            istart = ccnenx + nenx + sepnx;
            ncells = penx;
            cells = istart + (1 : ncells)';
            params.PositiveElectrode = NMC111Electrode('pe', G, cells);

            %% setup ccne
            istart = 1;
            ncells = ccnenx;
            cells = istart + (1 : ncells)';
            params.NegativeCurrentCollector = CurrentCollector('ccne', G, cells);

            %% setup ccpe
            istart = ccnenx + nenx + sepnx + penx;
            ncells = ccpenx;
            cells = istart + (1 : ncells)';
            params.PositiveCurrentCollector = CurrentCollector('ccpe', G, cells);

            %% setup sep
            istart = ccnenx + nenx;
            ncells = sepnx;
            cells = istart + (1 : ncells)';
            params.sep = celgard2500('sep', G, cells);
        
        end
        
        function coupTerm = setupNegativeCurrentCollectorBcCoupTerm(params)

            compnames = {'NegativeCurrentCollector'};
            coupTerm = couplingTerm('bc-NegativeCurrentCollector', compnames);
            coupTerm.couplingfaces = 1;
            coupTerm.couplingcells = 1;

        end

        function coupTerm = setupNegativeCurrentCollectorNegativeElectrodeCoupTerm(params)

            ccnenx = params.ccnenx;
            
            compnames = {'NegativeCurrentCollector', 'NegativeElectrode'};
            coupTerm = couplingTerm('NegativeCurrentCollector-NegativeElectrode', compnames);
            coupTerm.couplingfaces = [ccnenx + 1, 1];
            coupTerm.couplingcells = [ccnenx, 1];

        end
        
        function coupTerm = setupPositiveCurrentCollectorBcCoupTerm(params)

            ccpenx = params.ccpenx;

            compnames = {'PositiveCurrentCollector'};
            coupTerm = couplingTerm('bc-PositiveCurrentCollector', compnames);
            coupTerm.couplingfaces = ccpenx + 1;
            coupTerm.couplingcells = ccpenx;

        end

        function coupTerm = setupPositiveCurrentCollectorPositiveElectrodeCoupTerm(params)

            penx = params.penx;
            
            compnames = {'PositiveCurrentCollector', 'PositiveElectrode'};
            coupTerm = couplingTerm('PositiveCurrentCollector-PositiveElectrode', compnames);
            coupTerm.couplingfaces = [1, penx + 1];
            coupTerm.couplingcells = [1, penx];
            
        end

        function coupTerm = setupNegativeElectrodeElectrolyteCoupTerm(params)
            
            nenx   = params.nenx;

            compnames = {'NegativeElectrode', 'Electrolyte'};
            coupTerm = couplingTerm('NegativeElectrode-Electrolyte', compnames);
            cells1 = (1 : nenx)';
            cells2 = (1 : nenx)';
            coupTerm.couplingcells =  [cells1, cells2];
            coupTerm.couplingfaces = []; % no coupling throug faces. We set it as empty
            
        end
        
        function coupTerm = setupPositiveElectrodeElectrolyteCoupTerm(params)
            
            sepnx  = params.sepnx;
            nenx   = params.nenx;
            penx   = params.penx;
            
            compnames = {'PositiveElectrode', 'Electrolyte'};
            coupTerm = couplingTerm('PositiveElectrode-Electrolyte', compnames);
            cells1 = (1 : penx)';
            cells2 = nenx + sepnx + (1 : penx)';
            coupTerm.couplingcells = [cells1, cells2];
            coupTerm.couplingfaces = []; % no coupling between faces
            
        end

    end
    
end
