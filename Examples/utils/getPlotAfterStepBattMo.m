function fn = getPlotAfterStepBattMo(solver)
% Get a function that allows for dynamic plotting in `simulateScheduleAD`.
%
% SYNOPSIS:
%   fn = getPlotAfterStepBattMo(state0, model, schedule, 'plotWell', true);
%
%
% SEE ALSO:
%   `getPlotAfterStep`

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}


    h = figure();
    figureHandlers.voltageCurve = h;
    title('Cell Voltage / V');

    h = figure();
    figureHandlers.linearItPerNewtonIt = h;
    title('linear iterations/newton iterations');

    preconds = solver.LinearSolver.linearSolverSetup.preconditioners;
    nprecond = numel(preconds);
    for iprecond = 1 : nprecond
        h = figure();
        figureHandlers.preconds(iprecond) = h;
        titlestr = preconds(iprecond).name;
        titlestr = sprintf('preconditioner iteration/linear iteration - %s', titlestr);
        title(titlestr);        
    end

    drawnow;
    
    fn = @(model, states, reports, solver, schedule, simtime) afterStepFunction(model, states, reports, solver, schedule, simtime, figureHandlers);
    
end

function [model, states, reports, solver, ok] = afterStepFunction(model, states, reports, solver, schedule, simtime, figureHandlers)

    ok = true;
    
    setup = solver.LinearSolver.linearSolverSetup;
    
    fs = figureHandlers;
    
    ind = cellfun(@(x) ~isempty(x), reports);
    reports = reports(ind);

    ind = cellfun(@(x) not(isempty(x)), states); 
    states = states(ind);

    E = cellfun(@(x) x.Control.E, states); 
    I = cellfun(@(x) x.Control.I, states);
    time = cellfun(@(x) x.time, states); 


    setfig(fs.voltageCurve)
    plot(time, E, '*-');
    title('Cell Voltage / V')
    xlabel('time')

    
    its = getReportOutput(reports,'type','linearIterations');
    nits= getReportOutput(reports,'type','nonlinearIterations');

    setfig(fs.linearItPerNewtonIt)
    plot(its.time, its.total./nits.total, '*-');
    xlabel('time')
    title('linear iterations/newton iterations')

    switch setup.method

      case "grouped-gmres"

        % not implemented for the moment
        
      case "separate-variable-gmres"

        preconds = setup.preconditioners;
        npreconds = numel(preconds);
        pits = zeros(numel(its.time), npreconds);

        counter = 1;
        
        for icontrol = 1 : numel(reports)
            
            creport = reports{icontrol};

            for istep = 1 : numel(creport.StepReports)

                sreport = creport.StepReports{istep};
                
                for k = 1 : numel(sreport.NonlinearReport);
                    
                    nreport = sreport.NonlinearReport{k};
                    
                    for iprecond = 1 : npreconds

                        try
                            it = nreport.LinearSolver.precondReports{iprecond}.Iterations;
                            pits(counter, iprecond) = pits(counter, iprecond) + it;
                        end
                        
                    end
                end

                counter = counter + 1;
            end
        end
        
        for iprecond = 1 : npreconds
            
            setfig(fs.preconds(iprecond))
            plot(its.time, pits(:, iprecond)./its.total, '*-');
            xlabel('time');
            titlestr = preconds(iprecond).name;
            titlestr = sprintf('preconditioner iteration/linear iteration - %s', titlestr);
            title(titlestr);
            
        end
        

      otherwise
        
        error('setup.method not recognized');
    end

    drawnow;

end


function setfig(fignum)
    set(0, 'currentFigure', fignum);
end
