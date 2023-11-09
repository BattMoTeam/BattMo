classdef ExperimentalBattery < Battery
% 
% The battery model consists of 
%
% * an Electrolyte model given in :attr:`Electrolyte` property
% * a Negative Electrode Model given in :attr:`NegativeElectrode` property
% * a Positive Electrode Model given in :attr:`PositiveElectrode` property
% * a Thermal model given in :attr:`ThermalModel` property
%
    properties
        
        primaryVarNames
        funcCallList

    end
    
    methods
        
        function model = ExperimentalBattery(paramobj)

            model = model@Battery(paramobj);
            
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
        % Assembly of the governing equation
            
            opts = struct('ResOnly', false, 'iteration', 0, 'reverseMode', false); 
            opts = merge_options(opts, varargin{:});
            
            time = state0.time + dt;
            if(not(opts.ResOnly) && not(opts.reverseMode))
                state = model.initStateAD(state);
            elseif(opts.reverseMode)
               disp('No AD initatlization in equation old style')
               state0 = model.initStateAD(state0);
            else
                assert(opts.ResOnly);
            end


            %% We call the assembly equations ordered from the graph
            
            funcCallList = model.funcCallList;

            for ifunc = 1 : numel(funcCallList)
                eval(funcCallList{ifunc});
            end
            
            %% for now temperature and SOC are kept constant
            nc = model.G.cells.num;
            %state.SOC = model.SOC*ones(nc, 1);
            
            % Shorthands used in this function
            battery = model;
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            co      = 'Coating';
            am      = 'ActiveMaterial';
            cc      = 'CurrentCollector';
            elyte   = 'Electrolyte';
            am      = 'ActiveMaterial';
            itf     = 'Interface';
            sd      = "SolidDiffusion";
            thermal = 'ThermalModel';
            ctrl    = 'Control';
            
            %% We collect the governing equations
            % The governing equations are the mass and charge conservation equations for the electrolyte and the
            % electrodes and the solid diffusion model equations and the control equations. The equations are scaled to
            % a common value.

            massConsScaling = model.con.F;

            eqs = {};
            names = {};
            
            % Equation name : 'elyte_massCons';
            eqs{end + 1}   = state.(elyte).massCons*massConsScaling;
            names{end + 1} = 'elyte_massCons';

            % Equation name : 'elyte_chargeCons';
            eqs{end + 1}   = state.(elyte).chargeCons;
            names{end + 1} = 'elyte_chargeCons';
            
            % Equation name : 'ne_am_chargeCons';
            eqs{end + 1}   = state.(ne).(co).chargeCons;
            names{end + 1} = 'ne_co_chargeCons';

            % Equation name : 'pe_co_chargeCons';
            eqs{end + 1}   = state.(pe).(co).chargeCons;
            names{end + 1} = 'pe_co_chargeCons';
            
            switch model.(ne).(co).(am).diffusionModelType
                
              case 'simple'
                
                eqs{end + 1}   = state.(ne).(co).(am).(sd).massCons*massConsScaling;
                names{end + 1} = 'ne_co_am_sd_massCons';
                eqs{end + 1}   = state.(ne).(co).(am).(sd).solidDiffusionEq.*massConsScaling.*battery.(ne).(co).G.cells.volumes/dt;
                names{end + 1} = 'ne_co_am_sd_soliddiffeq';
                
              case 'full'
                
                % Equation name : 'ne_am_sd_massCons';
                n    = model.(ne).(co).(am).(itf).numberOfElectronsTransferred; % number of electron transfer (equal to 1 for Lithium)
                F    = model.con.F;
                vol  = model.(ne).(co).operators.pv;
                rp   = model.(ne).(co).(am).(sd).particleRadius;
                vsf  = model.(ne).(co).(am).(sd).volumetricSurfaceArea;
                surfp = 4*pi*rp^2;
                
                scalingcoef = (vsf*vol(1)*n*F)/surfp;
                
                eqs{end + 1}   = scalingcoef*state.(ne).(co).(am).(sd).massCons;
                names{end + 1} = 'ne_co_am_sd_massCons';
                eqs{end + 1}   = scalingcoef*state.(ne).(co).(am).(sd).solidDiffusionEq;
                names{end + 1} = 'ne_co_am_sd_soliddiffeq';

              otherwise
                
                error('diffusionModelType not recognized');
                
            end
            
            switch model.(pe).(co).(am).diffusionModelType
                
              case 'simple'

                eqs{end + 1}   = state.(pe).(co).(am).(sd).massCons*massConsScaling;
                names{end + 1} = 'pe_co_am_massCons';
                eqs{end + 1}   = state.(pe).(co).(am).(sd).solidDiffusionEq.*massConsScaling.*battery.(pe).(co).G.cells.volumes/dt;
                names{end + 1} = 'pe_co_am_sd_soliddiffeq';
                
              case 'full'
                
                % Equation name : 'pe_co_am_sd_massCons';
                n    = model.(pe).(co).(am).(itf).numberOfElectronsTransferred; % number of electron transfer (equal to 1 for Lithium)
                F    = model.con.F;
                vol  = model.(pe).(co).operators.pv;
                rp   = model.(pe).(co).(am).(sd).particleRadius;
                vsf  = model.(pe).(co).(am).(sd).volumetricSurfaceArea;
                surfp = 4*pi*rp^2;
                
                scalingcoef = (vsf*vol(1)*n*F)/surfp;

                eqs{end + 1} = scalingcoef*state.(pe).(co).(am).(sd).massCons;
                names{end + 1} = 'pe_co_am_sd_massCons';
                eqs{end + 1} = scalingcoef*state.(pe).(co).(am).(sd).solidDiffusionEq;
                names{end + 1} = 'pe_co_am_sd_soliddiffeq';
                
              otherwise
                
                error('diffusionModelType not recognized');
            end
            
            % Equation name : 'ne_cc_chargeCons';
            if model.(ne).include_current_collectors
                
                eqs{end + 1}   = state.(ne).(cc).chargeCons;
                names{end + 1} = 'ne_cc_chargeCons';
                
            end
            
            % Equation name : 'pe_cc_chargeCons';
            if model.(pe).include_current_collectors
                
                eqs{end + 1}   = state.(pe).(cc).chargeCons;
                names{end + 1} = 'pe_cc_chargeCons';
                
            end

            % Equation name : 'energyCons';
            if model.use_thermal
                
                eqs{end + 1}   = state.(thermal).energyCons;
                names{end + 1} = 'energyCons';
                
            end
            
            % Equation name : 'EIeq';
            eqs{end + 1}   = - state.(ctrl).EIequation;
            names{end + 1} = 'EIeq';
            
            % Equation name : 'controlEq'                                    
            eqs{end + 1}   = state.(ctrl).controlEquation;
            names{end + 1} = 'controlEq';

            neq = numel(eqs);
            types = repmat({'cell'}, 1, neq);
            primaryVars = model.getPrimaryVariableNames();

            %% Setup LinearizedProblem that can be processed by MRST Newton API
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
            
        end

        function primaryvarnames = getPrimaryVariableNames(model)

            primaryvarnames = model.primaryVarNames;
            
        end
        

        function model = validateModel(model, varargin)

            if isempty(model.computationalGraph)
                model = model.setupComputationalGraph();
                cgt = model.computationalGraph;
                model.primaryVarNames = cgt.getPrimaryVariableNames();
                model.funcCallList = cgt.getOrderedFunctionCallList();
            end
            
        end
        

    end

    
end



%{
Copyright 2021-2023 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
