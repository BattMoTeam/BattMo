classdef ProtonicMembraneCellWithGasSupply < BaseModel
    
    properties
        
        % Temperature
        T
        % Structure with physical constants
        constants

        Cell
        GasSupply

        couplingTerms
        couplingnames

        funcCallList
        primaryVarNames
        equationVarNames
        equationNames
        equationTypes

        scalings
        
    end
    
    methods
        
        function model = ProtonicMembraneCellWithGasSupply(paramobj)

            model = model@BaseModel();

            fdnames = {'T'       , ...
                       'couplingTerms'};
            
            model = dispatchParams(model, paramobj, fdnames);

            model.Cell      = ProtonicMembraneCell(paramobj.Cell);
            model.GasSupply = ProtonicMembraneGasSupply(paramobj.GasSupply);
            
            model.constants = PhysicalConstants();
            
        end
        
        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@BaseModel(model);

            fn =  @ProtonicMembraneCellWithGasSupply.updateBcAnodePressures;
            inputvarnames = {VarName({'GasSupply', 'GasSupplyBc'}, 'pressures', 2)};
            outputvarname = VarName({'Cell', 'Anode'}, 'pressures', 2);
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});
            
        end

        function [problem, state] = getEquations(model, state0, state,dt, drivingForces, varargin)
            
            sd  = 'SolidDiffusion';
            itf = 'Interface';
            
            time = state0.time + dt;
            state = model.initStateAD(state);
            
            %% We call the assembly equations ordered from the graph

            funcCallList = model.funcCallList;

            for ifunc = 1 : numel(funcCallList)
                eval(funcCallList{ifunc});
            end

            state = model.applyScaling(state);
            
            for ieq = 1 : numel(model.equationVarNames)
                eqs{ieq} = model.getProp(state, model.equationVarNames{ieq});
            end
            
            names       = model.equationNames;
            types       = model.equationTypes;
            primaryVars = model.primaryVarNames;
            
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

        end

        function primaryvarnames = getPrimaryVariableNames(model)

            primaryvarnames = model.primaryVarNames;

        end
        
    end
    
end
