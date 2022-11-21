function stepreports = cleanupReport(report, varargin)

    opt = struct('saveToJson', false, ...
                 'jsonFileName', []);

    opt = merge_options(opt, varargin{:});

    if isfield(report, 'ControlstepReports')
        report = report.ControlstepReports;
    end

    for ireport = 1 : numel(report)
         report{ireport} = cleanupControlTimeReport(report{ireport});
    end

    if opt.saveToJson
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

    for istepreport = 1 : numel(ctreport.StepReports)
        ctreport.StepReports{istepreport} = cleanupStepReport(ctreport.StepReports{istepreport});
        ctreport.StepReports{istepreport}.index = istepreport;
    end
    
end


function stepreport = cleanupStepReport(stepreport)

    
    stepreport = rmfield(stepreport, 'Timestep');
    stepreport = rmfield(stepreport, 'LocalTime');

    stepreport.NonlinearReport = cleanupNonlinearReport(stepreport.NonlinearReport);
    
end

function nonlinearReport = cleanupNonlinearReport(nonlinearReport)

    for inr = 1 : numel(nonlinearReport)

        nonlinearReport{inr} = rmfield(nonlinearReport{inr}, 'UpdateState');
        nonlinearReport{inr} = rmfield(nonlinearReport{inr}, 'Failure');
        nonlinearReport{inr} = rmfield(nonlinearReport{inr}, 'FailureMsg');
        nonlinearReport{inr} = rmfield(nonlinearReport{inr}, 'FinalUpdate');
        nonlinearReport{inr} = rmfield(nonlinearReport{inr}, 'Residuals');
        nonlinearReport{inr} = rmfield(nonlinearReport{inr}, 'AssemblyTime');
        nonlinearReport{inr} = rmfield(nonlinearReport{inr}, 'StabilizeReport');
        nonlinearReport{inr} = rmfield(nonlinearReport{inr}, 'ResidualsConverged');

        nonlinearReport{inr}.index = inr;
        nonlinearReport{inr}.LinearSolver = cleanupLinearSolverReport(nonlinearReport{inr}.LinearSolver);
    end
    
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
    
end

    
