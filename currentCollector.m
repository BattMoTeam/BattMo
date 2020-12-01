classdef currentCollector < FvModel
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        % Physical constants
        con = physicalConstants();
        
        % Design properties
        t       % Thickness,        [m]
        
        % Material properties
        am % Active material object
        
        % State properties
        E % Electric potential,   [V]
        j % Current density,      [A/m2]
        T % Temperature
        
        % Mesh properties
        X
        Xb
        N
        dombin
                
    end
    
    method
        function obj = currentCollector(T)
            
            obj.am  = currentCollectorAM(T);
            
            obj.am.eps = 1;
            
            obj.eps = 1;
            
            obj.t = 88e-6;
            
            obj.T = T;
            
            obj.thermodynamics()
            obj.E = obj.am.OCP;
        end
        
        function thermodynamics(obj)
            %THERMODYNAMICS Summary of this method goes here
            %   Detailed explanation goes here
           
            
            % Electrochemical Thermodynamics
            obj.am.equilibrium();
            
            
        end
        

        
    end
end

