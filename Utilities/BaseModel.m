classdef BaseModel < PhysicalModel
    
    methods

        function model = BaseModel()
            model = model@PhysicalModel([]);
        end
        
        function state = setProp(model, state, names, val)
            if iscell(names) & (numel(names) > 1)
                name = names{1};
                names = names(2 : end);
                state.(name) = model.setProp(state.(name), names, val);
            elseif iscell(names) & (numel(names) == 1)
                name = names{1};
                if isnumeric(name)
                    state{name} = val;
                else
                    state.(name) = val;
                end
            else
                error('format not recognized');
            end
        end
        
        
        function state = setNewProp(model, state, names, val)
            if iscell(names) & (numel(names) > 1)
                name = names{1};
                names = names(2 : end);
                if(not(isfield(state,name)))
                    state.(name)=struct();
                end
                state.(name) = model.setNewProp(state.(name), names, ...
                                                               val);
            elseif iscell(names) & (numel(names) == 1)
                name = names{1};
                if isnumeric(name)
                    if(not(iscell(state)))
                        state = {};
                    end
                    state{name} = val;
                else
                    state.(name) = val;
                end
            else
                error('format not recognized');
            end
        end

        
        function var = getProp(model, state, names)
            if iscell(names) && (numel(names) > 1)
                name = names{1};
                names = names(2 : end);
                var = model.getProp(state.(name), names);
            elseif iscell(names) & (numel(names) == 1)
                name = names{1};
                if isnumeric(name)
                    var = state{name};
                else
                    var = state.(name);
                end
            else
                error('format not recognized');
            end
        end
        
        function submod = getSubmodel(model, names)
            submod = model.(names{1});
            for i=2:numel(names)
                submod = submod.(names{i});
            end
        end


        function [state, report] = updateState(model, state, problem, dx, drivingForces)

            p = model.getPrimaryVariables();

            for i = 1 : numel(dx)
                val = model.getProp(state, p{i});
                val = val + dx{i};
                state = model.setProp(state, p{i}, val);
            end
           
            report = [];
        end
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)

            [state, report] = updateAfterConvergence@PhysicalModel(model, state0, state, dt, drivingForces);
            p = model.getPrimaryVariables();
            cleanState = [];
            for ind = 1 : numel(p)
                cleanState = model.copyProp(cleanState, state, p{ind});
            end
            cleanState.time = state.time;
            
            state = cleanState;
            report = [];
            
        end
        
        function state = copyProp(model, state, refState, names)
            
            if iscell(names) & (numel(names) > 1)
                name = names{1};
                names = names(2 : end);
                if isfield(state, name)
                    state.(name) = model.copyProp(state.(name), refState.(name), names);
                else
                    state.(name) = model.copyProp([], refState.(name), names);
                end
            elseif iscell(names) & (numel(names) == 1)
                name = names{1};
                if isnumeric(name)
                    state{name} = refState{name};
                else
                    state.(name) = refState.(name);
                end
            else
                error('format not recognized');
            end
        end

        
        function scale = getScales()
            scale = [];
            error();
        end
        
        function [state, report] = updateStateNew(model, state, problem, ...
                                                  dx, drivingForces)
            scales = model.getScales();
            for i = 1:numel(problem.primaryVariables)
                p = problem.primaryVariables{i};
                % Update the state
                scale=model.getProp(scales,p);
                if(isempty(scale))
                     state = model.updateStateFromIncrement(state, dx{i}, ...
                                                            problem, p);
                else
                    state = model.updateStateFromIncrement(state, dx{i}, ...
                                                           problem, p, ...
                                                           scale.relchangemax, ...
                                                           scale.abschangemax);
                    val = model.getProp(state, p);
                    val = max(val,scale.min);
                    val = min(val,scale.max);
                    state = model.setProp(state, p, val);
                end               
            end
            report = []
        end
        
        function state = reduceState(model, state, removeContainers)
            state = value(state);
        end
        
    end
    
end
