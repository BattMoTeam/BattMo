classdef SeaWaterBatteryGeneratorP2D < SeaWaterBatteryGenerator
% setup 1D grid 
    properties
        
        andenx = 10; % number of cells in anode
        ctdenx = 10; % number of cells in cathode 
        sepnx  = 10; % number of cells in electrolyte between the cathode and anode 
        fac    = 1;
        
    end
    
    methods
        
        function gen = SeaWaterBatteryGeneratorP2D()
          gen = gen@SeaWaterBatteryGenerator();  
        end
            
        function inputparams = updateBatteryInputParams(gen, inputparams)
            inputparams = gen.setupBatteryInputParams(inputparams, []);
        end
        
        function [inputparams, gen] = setupGrid(gen, inputparams, ~)
        % inputparams is instance of BatteryInputParams
        % setup inputparams.G
        
            andenx = gen.andenx;
            ctdenx = gen.ctdenx;
            sepnx  = gen.sepnx;

            nxs = [andenx; sepnx; ctdenx];

            % default value taken from original SeaWater code 
            xlength = 1e-6*[5000; 1000; 100];

            x = xlength./nxs;
            x = rldecode(x, nxs);
            x = [0; cumsum(x)];

            G = tensorGrid(x);
            G = computeGeometry(G); 

            inputparams.G = G;
            gen.G = G;
            
        end

        function gen = applyResolutionFactors(gen)
            
            fac = gen.fac;
            
            gen.andenx = gen.andenx*fac;
            gen.ctdenx = gen.ctdenx*fac;
            gen.sepnx  = gen.sepnx*fac;
            
        end
            
        function inputparams = setupElectrolyte(gen, inputparams, params)
            
            % In this case we setup the electrolyte as a subgrid of the background, even if it the two in fact
            % coincides. It is unecessary but we do it to keep the approach generic.
            params.cellind = (1 : (gen.andenx + gen.sepnx + gen.ctdenx))';
            inputparams = setupElectrolyte@SeaWaterBatteryGenerator(gen, inputparams, params);
        end
        
        function inputparams = setupElectrodes(gen, inputparams, params)

            ande  = 'Anode';
            ctde  = 'Cathode';
            
            andenx = gen.andenx; 
            ctdenx = gen.ctdenx; 
            sepnx  = gen.sepnx; 
            
            %% parameters for anode
            params.(ande).cellind = (1 : andenx)';
            params.(ande).bcfaces = 1;
            params.(ande).bccells = 1;
            
            %% parameters for cathode
            params.(ctde).cellind = (andenx + sepnx) + (1 : ctdenx)';
            params.(ctde).bcfaces = ctdenx + 1;
            params.(ctde).bccells = ctdenx;
            
            inputparams = setupElectrodes@SeaWaterBatteryGenerator(gen, inputparams, params);
            

        end            I
                
    end
    
end


