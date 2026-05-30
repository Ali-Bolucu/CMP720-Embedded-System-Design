%% RunStatisticalBenchmark.m
% Runs multiple independent trials to provide statistical reliability
% (Mean, Std, Median) for the Spectral Seeding vs Pure Random comparison.

function RunStatisticalBenchmark()
    clc;
    fprintf('==============================================================\n');
    fprintf('      STATISTICAL BENCHMARK (MULTIPLE RUNS) STARTING\n');
    fprintf('==============================================================\n');
    
    addpath(fileparts(mfilename('fullpath')));
    
    % --- CONFIGURATION ---
    numRuns = 30;  % Standard in academia for statistical significance
    Npop = 100;
    MaxFE = 10000;
    
    fprintf('Configuration:\n');
    fprintf('  Number of Independent Runs : %d\n', numRuns);
    fprintf('  Population Size            : %d\n', Npop);
    fprintf('  Max Evaluations            : %d\n\n', MaxFE);
    
    % Initialize MegaGraph
    MegaGraph = ApplicationToMegaGraph();
    
    % Storage for Results
    results_random = struct('hv', zeros(numRuns,1), 'igd', zeros(numRuns,1), 'time', zeros(numRuns,1), 'bestLat', zeros(numRuns,1), 'bestTh', zeros(numRuns,1));
    results_spectral = struct('hv', zeros(numRuns,1), 'igd', zeros(numRuns,1), 'time', zeros(numRuns,1), 'bestLat', zeros(numRuns,1), 'bestTh', zeros(numRuns,1));
    
    % Storage for all fronts to calculate global IGD
    all_fronts_random = cell(numRuns, 1);
    all_fronts_spectral = cell(numRuns, 1);
    
    global USE_SPECTRAL_SEEDING
    
    for runIdx = 1:numRuns
        fprintf('--- RUN %d / %d ---\n', runIdx, numRuns);
        
        % ---------------------------------------------------------
        % 1. PURE RANDOM APPROACH
        % ---------------------------------------------------------
        USE_SPECTRAL_SEEDING = false;
        
        t1 = tic;
        [~, Obj_random, ~] = platemo('N', Npop, 'maxFE', MaxFE, 'problem', @NocMappingProblem, 'algorithm', @NSGAII, 'parameter', MegaGraph, 'save', 0);
        results_random.time(runIdx) = toc(t1);
        
        % Extract non-dominated and calculate stats
        FrontObj_r = GetNonDominated(Obj_random);
        all_fronts_random{runIdx} = FrontObj_r;
        results_random.bestLat(runIdx) = min(FrontObj_r(:, 1));
        results_random.bestTh(runIdx)  = min(FrontObj_r(:, 2));
        
        % ---------------------------------------------------------
        % 2. SPECTRAL SEEDING APPROACH
        % ---------------------------------------------------------
        USE_SPECTRAL_SEEDING = true;
        
        t2 = tic;
        [~, Obj_spectral, ~] = platemo('N', Npop, 'maxFE', MaxFE, 'problem', @NocMappingProblem, 'algorithm', @NSGAII, 'parameter', MegaGraph, 'save', 0);
        results_spectral.time(runIdx) = toc(t2);
        
        % Extract non-dominated and calculate stats
        FrontObj_s = GetNonDominated(Obj_spectral);
        all_fronts_spectral{runIdx} = FrontObj_s;
        results_spectral.bestLat(runIdx) = min(FrontObj_s(:, 1));
        results_spectral.bestTh(runIdx)  = min(FrontObj_s(:, 2));
        
        % ---------------------------------------------------------
        % HV Calculation for this run
        % ---------------------------------------------------------
        % Calculate dynamic ref point for this run to compare them fairly
        refPoint = max([FrontObj_s; FrontObj_r]) .* 1.1;
        results_random.hv(runIdx) = Calculate2DHV(FrontObj_r, refPoint);
        results_spectral.hv(runIdx) = Calculate2DHV(FrontObj_s, refPoint);
        
        fprintf('  Random   -> HV: %.4f | Time: %.2fs\n', results_random.hv(runIdx), results_random.time(runIdx));
        fprintf('  Spectral -> HV: %.4f | Time: %.2fs\n', results_spectral.hv(runIdx), results_spectral.time(runIdx));
    end
    
    % ---------------------------------------------------------
    % GLOBAL IGD CALCULATION (Post-Processing)
    % ---------------------------------------------------------
    fprintf('\nPooling all fronts to create Global Reference Front for IGD...\n');
    pooled_fronts = [];
    for i = 1:numRuns
        pooled_fronts = [pooled_fronts; all_fronts_random{i}; all_fronts_spectral{i}];
    end
    globalRefFront = GetNonDominated(pooled_fronts);
    
    for i = 1:numRuns
        results_random.igd(i) = CalculateIGD(all_fronts_random{i}, globalRefFront);
        results_spectral.igd(i) = CalculateIGD(all_fronts_spectral{i}, globalRefFront);
    end
    
    % --- PRINT STATISTICAL SUMMARY ---
    fprintf('\n========================================================================================\n');
    fprintf('                               STATISTICAL SUMMARY (%d RUNS)\n', numRuns);
    fprintf('========================================================================================\n');
    fprintf('%-20s | %-12s | %-12s | %-12s | %-12s\n', 'Metric', 'Approach', 'Mean', 'Median', 'Std Dev');
    fprintf('----------------------------------------------------------------------------------------\n');
    
    PrintStatRow('Hypervolume', 'Pure Random', results_random.hv);
    PrintStatRow('', 'Spectral', results_spectral.hv);
    fprintf('----------------------------------------------------------------------------------------\n');
    PrintStatRow('IGD', 'Pure Random', results_random.igd);
    PrintStatRow('', 'Spectral', results_spectral.igd);
    fprintf('----------------------------------------------------------------------------------------\n');
    PrintStatRow('Best Latency', 'Pure Random', results_random.bestLat);
    PrintStatRow('', 'Spectral', results_spectral.bestLat);
    fprintf('----------------------------------------------------------------------------------------\n');
    PrintStatRow('Best Thermal', 'Pure Random', results_random.bestTh);
    PrintStatRow('', 'Spectral', results_spectral.bestTh);
    fprintf('----------------------------------------------------------------------------------------\n');
    PrintStatRow('Exec Time (s)', 'Pure Random', results_random.time);
    PrintStatRow('', 'Spectral', results_spectral.time);
    fprintf('========================================================================================\n');
    
    % Save to a MAT file for plotting later
    save('statistical_results.mat', 'results_random', 'results_spectral');
    fprintf('Detailed results saved to "statistical_results.mat".\n');
end

% --- Helper Functions ---

function PrintStatRow(metricName, approachName, dataArray)
    fprintf('%-20s | %-12s | %-12.4f | %-12.4f | %-12.4f\n', ...
            metricName, approachName, mean(dataArray), median(dataArray), std(dataArray));
end

function NonDom = GetNonDominated(Obj)
    % Simple non-dominated sorting for 2D objectives (Minimization)
    N = size(Obj, 1);
    isDominated = false(N, 1);
    for i = 1:N
        for j = 1:N
            if i ~= j
                % j dominates i if j is <= in all objectives and < in at least one
                if all(Obj(j,:) <= Obj(i,:)) && any(Obj(j,:) < Obj(i,:))
                    isDominated(i) = true;
                    break;
                end
            end
        end
    end
    NonDom = Obj(~isDominated, :);
end

function hv = Calculate2DHV(Front, RefPoint)
    % Very simple 2D Hypervolume calculation (minimization)
    % Sort front by first objective
    [~, idx] = sort(Front(:,1));
    sortedFront = Front(idx, :);
    
    hv = 0;
    prev_x = RefPoint(1);
    for i = 1:size(sortedFront, 1)
        cur_x = sortedFront(i, 1);
        cur_y = sortedFront(i, 2);
        
        width = prev_x - cur_x;
        height = RefPoint(2) - cur_y;
        
        if width > 0 && height > 0
            hv = hv + (width * height);
        end
        prev_x = cur_x;
    end
end

function igd = CalculateIGD(Front, RefFront)
    % Calculates the Inverted Generational Distance (IGD)
    if isempty(Front) || isempty(RefFront)
        igd = Inf;
        return;
    end
    distances = pdist2(RefFront, Front);
    minDistances = min(distances, [], 2);
    igd = mean(minDistances);
end
