classdef BatteryInputParams
    
    properties
        
        % Temperature and SOC
        T
        SOC

        % Input current
        J
        
        % Voltage cut
        Ucut
        
        % Grid of the whole model (parent grid)
        G

        %% Electrolyte model
        elyte
        
        %% Negative electrode model
        ne
        
        %% Positive electrode model
        pe
        
        %% Separator model
        sep
        
        %% Current collector for negative electrode model
        ccne
        
        %% Current collector for positive electrode model 
        ccpe
        
        %% Coupling terms (describe the topological structure of the coupling)

        % electrolyte - negative electrode (location where chemical reactions take place)
        coupTermNeElyte
        
        % electrolyte - positive electrode (location where chemical reactions take place)
        coupTermPeElyte 

        % current collector and electrode (current continuity is imposed at the faces)
        coupTermCcneNe  
        coupTermCcpePe  

        % Boundary conditions for the current collectors
        coupTermCcneBc  
        coupTermCcpeBc  
        
    end
    
    methods
        
        function params = BatteryInputParams(params)

            params = params.setupVariousParams();
            params = params.setupSubModels();
            params = params.setupCouplingTerms();
            
        end
        
        function params = setupVariousParams(param)
        % In this function, we define :  SOC, T,  J and  Ucut
            error('virtual function')
        end
        
        function params = setupSubModels(params)
        % In this function, we setup: G and the models (elyte, sep, ne, pe, ccne, ccpe)
            error('Virtual Function');
        end
        
        function params = setupCouplingTerms(params)
            params.coupTermNeElyte = setupNeElyteCoupTerm(params);
            params.coupTermPeElyte = setupPeElyteCoupTerm(params);
            params.coupTermCcneNe  = setupCcneNeCoupTerm(params);
            params.coupTermCcpePe  = setupCcpePeCoupTerm(params);
            params.coupTermCcneBc  = setupCcneBcCoupTerm(params);
            params.coupTermCcpeBc  = setupCcpeBcCoupTerm(params);
        end

        function coupTerm = setupCcneBcCoupTerm(params)
        % In this function, we setup coupling corresponding coupling 
            error('Virtual Function');
        end

        function coupTerm = setupCcneNeCoupTerm(params)
        % In this function, we setup coupling corresponding coupling 
            error('Virtual Function');
        end
        
        function coupTerm = setupCcpeBcCoupTerm(params)
        % In this function, we setup coupling corresponding coupling 
            error('Virtual Function');
        end

        function coupTerm = setupCcpePeCoupTerm(params)
        % In this function, we setup coupling corresponding coupling 
            error('Virtual Function');            
        end

        function coupTerm = setupNeElyteCoupTerm(params)
        % In this function, we setup coupling corresponding coupling 
            error('Virtual Function');
        end

        function coupTerm = setupPeElyteCoupTerm(params)
        % In this function, we setup coupling corresponding coupling 
            error('Virtual Function');
        end

    end
    
end
