%% ROBUSTNESSTEST for FittingEIS  

num_runs = 5;

%% PERTURBATION DISTANCE FACTOR 
% Controls the maximum relative deviation from the reference parameters.
% 0.05 = Tight search within +/- 5% of the reference
% 0.20 = Medium search within +/- 20% of the reference
% 0.60 = Wide search within +/- 60% of the reference
perturbation_factor = 0.20; 

% Load the full experimental dataset
eisdata = load_chen_data();

%%
% We define the reference parameters (the center point for our test)

params0 = [0.05052; 1.12673; 59119.9; 0.03155; 11054.0];  
feis = FittingEIS(eisdata, params0);
[best_params_ref, fitting_error] = feis.run();

% Extract hard boundary constraints from the class properties
pmin = feis.scales(1:5);
pmax = feis.scales(6:10);

% Pre-allocate logging arrays to store optimization history results
all_best_params = zeros(num_runs, 5);
all_errors = zeros(num_runs, 1);

fprintf('\n==================================================\n');
fprintf('    STARTING EIS OPTIMIZATION ROBUSTNESS TEST     \n');
fprintf('==================================================\n');
fprintf('Number of runs: %d | Perturbation Radius: +/-%.1f%%\n', num_runs, perturbation_factor * 100);

% --- MAIN ROBUSTNESS LOOP ---
for i = 1:num_runs
    
    % Generate a random variation vector between -1 and +1
    rand_deviation = 2 * rand(5, 1) - 1;
    
    % Scale the reference parameters by the perturbation factor
    params0_perturbed = best_params_ref(:) .* (1 + perturbation_factor * rand_deviation);
    
    % Enforce boundary box constraints so guesses don't breach pmin/pmax
    params0_perturbed = max(pmin(:), min(pmax(:), params0_perturbed));
    
    fprintf('\n-----------------------------------\n');
    fprintf('Run %d/%d - Perturbed Guess:\n', i, num_runs);
    disp(params0_perturbed'); % Displayed horizontally for cleaner logs
    
    try
        % 2. Instantiate the FittingEIS class with the new randomized guess
        feis = FittingEIS(eisdata, params0_perturbed);
        
        % 3. Run the objective optimization process 
        
        % 4. Log the output metrics for the current evaluation run
        all_best_params(i, :) = best_params(:)';
        all_errors(i) = fitting_error;
        
        % Display run metrics inside the console
        fprintf('-> Success Run %d! Final Error (SSR): %.4e\n', i, fitting_error);
        disp('Optimized parameters found:');
        disp(best_params');
        
    catch ME
        % If a specific initial guess causes numerical instability, skip it safely
        fprintf('-> Run %d failed: %s\n', i, ME.message);
        all_best_params(i, :) = NaN;
        all_errors(i) = NaN;
    end
    
end

%% --- GLOBAL PERFORMANCE SUMMARY ---

fprintf('\n==================================================\n');
fprintf('               ROBUSTNESS SUMMARY                 \n');
fprintf('==================================================\n');
for i = 1:num_runs
    if any(isnan(all_best_params(i, :)))
        fprintf('Run %d : FAILED\n', i);
    else
        fprintf('Run %d : SSR Error = %.4e | R0 = %.2e | R1 = %.2e | C1 = %.2e\n', ...
            i, all_errors(i), all_best_params(i,1), all_best_params(i,2), all_best_params(i,3));
    end
end
