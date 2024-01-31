function stepreports = cleanupReport(report, varargin)

    opt = struct('saveToJson', false, ...
                 'jsonFileName', []);

    opt = merge_options(opt, varargin{:});

    if isfield(report, 'ControlstepReports')
        report = report.ControlstepReports;
    end

    for ireport = 1 : numel(report)
        report{ireport}.index = ireport;
        report{ireport} = cleanupControlTimeReport(report{ireport});
    end

    if opt.saveToJson | ~isempty(opt.jsonFileName)
        jsonFileName = opt.jsonFileName;
        if isempty(jsonFileName)
            jsonFileName = 'temp.json';
        end
        fid = fopen(jsonFileName, 'w');
        jsonstr = jsonencode(report);
        fprintf(fid, jsonstr);
        fclose(fid);
    end

end

function ctreport = cleanupControlTimeReport(ctreport)

    fds = {'Iterations', ...
           'EarlyStop' , ...
           'WallTime'  , ...
           'MinistepCuttingCount'};
    for ifd = 1 : numel(fds)
        ctreport = rmfield(ctreport, fds{ifd});
    end

    ctreport = setfieldfirst(ctreport, 'index');
    
    for istepreport = 1 : numel(ctreport.StepReports)
        ctreport.StepReports{istepreport}.index = istepreport;
        ctreport.StepReports{istepreport} = cleanupStepReport(ctreport.StepReports{istepreport});
    end

end


function stepreport = cleanupStepReport(stepreport)

    
    stepreport = rmfield(stepreport, 'Timestep');
    stepreport = rmfield(stepreport, 'LocalTime');

    stepreport.nonlinearIterations = numel(stepreport.NonlinearReport);
    
    stepreport = setfieldfirst(stepreport, 'nonlinearIterations');
    stepreport = setfieldfirst(stepreport, 'index');
    
    stepreport.NonlinearReport = cleanupNonlinearReport(stepreport.NonlinearReport);
    
end

function nonlinearReport = cleanupNonlinearReport(nonlinearReport)

    for inr = 1 : numel(nonlinearReport)
        nonlinearReport{inr}.index = inr;
        nonlinearReport{inr} = cleanupNonlinearReportItem(nonlinearReport{inr});
    end
    
end

function nonlinearReportItem = cleanupNonlinearReportItem(nonlinearReportItem)

    fds = {'UpdateState'    , ...
           'Failure'        , ...
           'FailureMsg'     , ...
           'FinalUpdate'    , ...
           'Residuals'      , ...
           'AssemblyTime'   , ...
           'StabilizeReport', ...
           'ResidualsConverged'};

    for ifd = 1 : numel(fds)
        nonlinearReportItem = rmfield(nonlinearReportItem, fds{ifd});
    end

    nonlinearReportItem = setfieldfirst(nonlinearReportItem, 'index');
    
    nonlinearReportItem.LinearSolver = cleanupLinearSolverReport(nonlinearReportItem.LinearSolver);
    
end
    

function lreport = cleanupLinearSolverReport(lreport)
    fdnames = {'Residual',
               'SolverTime',
               'failureMsg',
               'LinearSolutionTime',
               'PreparationTime',
               'PostProcessTime',
               'Converged'};

    for ifd = 1 : numel(fdnames)
        if isfield(lreport, fdnames{ifd})
            lreport = rmfield(lreport, fdnames{ifd});
        end
    end

    lreport = setfieldfirst(lreport, 'index');
    
end

    
function s = setfieldfirst(s, fd)
    fds = fieldnames(s);
    ind = ismember(fds, fd);
    ind = [find(ind); find(~ind)];
    s = orderfields(s, ind);
end



%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
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
