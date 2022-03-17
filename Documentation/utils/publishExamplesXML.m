function publishExamplesXML(varargin)
%
% Example: publishExamplesXML('exampleNames', {'runBattery1D'})
% 
% The scripts ('runBattery1D' in example above) should be in path
%
% To convert the examples from xml to rst format (as needed by documentation), it should be enough to run 
%    python buildPublishedExamples
% in terminal (script buildPublishedExamples.py is located in the same directory as publishExamplesXML.m)
    opt = struct('replace'     , false, ...
                 'extension'   , 'xml', ...
                 'catchError'  , false, ...
                 'exampleNames', []);
    opt = merge_options(opt, varargin{:});
    
    exampleNames = opt.exampleNames;

    if isempty(exampleNames)
        error('Please, provide list of examples you want to publish from the Examples directory');
    end
        
    BattMoExampleDir = mfilename('fullpath'); 
    BattMoExampleDir = fileparts(BattMoExampleDir);
    BattMoExampleDir = fullfile(BattMoExampleDir, '..', '..', 'Examples');

    outputDir = fullfile(BattMoExampleDir, '..', 'Documentation', 'publishedExamples');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end 
    
    mrstVerbose off
    
    count = 0;
    total = 0;
    
    for exNo = 1:numel(exampleNames)
        exampleName = exampleNames{exNo};
        filename = fullfile(BattMoExampleDir, exampleName);
        
        isOk = isExamplePublished(filename, opt.extension);

        if ~isOk || opt.replace
            fprintf('Publishing example %s...', exampleName);
            close all;
            pstatus = pause('query');
            pause('off')
            try
                publish_opt = struct('format'        , opt.extension , ...
                                     'maxOutputLines', 8             , ...
                                     'catchError'    , opt.catchError, ....
                                     'outputDir'     , outputDir);
                run_publish(filename, publish_opt)
                fprintf(' Done.\n');
                count = count + 1;
            catch e
                fprintf('\n *** Error in publish example: %s\n', e.message);
            end
            pause(pstatus);
            close all; 
        else
            fprintf('Example %s already published, skipping...\n', exampleName);
        end
        
        total = total + 1;
        
    end
    fprintf('All done! Published %d of %d total examples\n', count, total);
    
end

function [pth, name, ext] = getPublishedPath(filename)
    
    [expth, filename, ext] = fileparts(filename);
    pth = fullfile(expth, 'html');
    name = filename;
    
end

function isPublished = isExamplePublished(filename, ext)
    
    [pth, name, ~] = getPublishedPath(filename);
    isPublished = exist(fullfile(pth, [name, '.', ext]), 'file') > 0;
    
end

function run_publish(filename, opt)
    
    publish(filename, opt);
    
end

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
