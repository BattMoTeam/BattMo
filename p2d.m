classdef p2d < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        % Physical constants
        con = physicalConstants();
        
        % Component objects
        ne
        pe
        sep
        elyte
        
        % Finite volume structure
        fv
        
        
    end
    
    methods
        function obj = p2d(inputArg1,inputArg2)
            %P2D Construct an instance of the P2D model class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function buildMesh(obj)
            
            
        end
        
        function ic(obj)
            %IC Sets the initial conditions for the P2D simulation
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
        
        function dynamicPreprocess(obj)
            
            
        end
        
        function readState(obj)
            
            
        end
        
        
        function buildSOE(obj)
            
            
        end
        
        
        
    end
end

