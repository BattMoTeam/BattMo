classdef fv1d
    %FV1D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        % Domain identification
        domnames    % Domain names
        domid       % Domain ID
        dombin      % Domain binary vectors
        dom         % Domain data structures
        
        % Mesh properties
        dXvec       % Nominal cell size vector
        Lvec        % Nominal domain length vector
        Nvec        % Nominal domain cell quantity vector
        N           % Total cell quantity
        dX          % Cell size vector
        Xb          % Cell boundary locations
        X           % Cell node locations
        
        % Time discretization properties
        ti          % Initial time, s
        tf          % Final time, s
        dt          % Time step, s
        tUp         % Ramp up time, s
        tSpan       % Time span vector, s
        
        % State properties
        y
        y0
        yp
        yp0
        
        % State data slots
        s1, s2, s3, s4, s5
        s6, s7, s8, s9, s10
        s11, s12, s13, s14, s15
        s16, s17, s18, s19, s20
        
    end
    
    methods
        function obj = fv1d(names, sizes, lengths)
            %FV1D Construct an instance of this class
            %   Detailed explanation goes here
            %% Finite volume mesh properties
            % Define finite volume domains, cell sizes, component lengths
            % and number of cells in the mesh
            obj.domnames = names;
            obj.dXvec    = sizes;
            obj.Lvec     = lengths;
           
            obj.Nvec    = round(obj.Lvec ./ obj.dXvec);
            obj.N = sum(obj.Nvec);
                        
            % Build FV size vector
            obj.dX = [];
            for i = 1:length(obj.Lvec)
                obj.dX = vertcat(obj.dX, obj.dXvec(i) .* ones(obj.Nvec(i),1));
            end
            
            % Build center and boundary node vector
            obj.Xb = zeros(obj.N+1,1);
            obj.X = zeros(obj.N,1);
            for i = 2:length(obj.dX)+1
                obj.Xb(i) = obj.Xb(i-1) + obj.dX(i-1);
                obj.X(i-1) = obj.Xb(i-1) + 0.5*obj.dX(i-1);
            end
            
            % Define mesh slots and binary vectors
            for i = 1:length(obj.Lvec)
                if i == 1
                    obj.domid{i} = 1:obj.Nvec(i);
                else
                    obj.domid{i} = (1+obj.domid{i-1}(end)):(obj.Nvec(i)+obj.domid{i-1}(end));
                end
                
                obj.dombin{i} = zeros(obj.N,1);
                obj.dombin{i}(obj.domid{i}(:)) = 1; 
            end        
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

