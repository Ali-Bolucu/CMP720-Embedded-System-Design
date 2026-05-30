%% RunAlgorithmComparison.m
% Tests the Spectral Seeding initialization against multiple MOEAs
% to prove that the seeding method is universally beneficial, not just
% limited to NSGA-II.

function RunAlgorithmComparison()
    clc;
    fprintf('=========================================================================\n');
    fprintf('           ALGORITHMIC COMPARISON (NSGA-II vs NSGA-III vs SPEA2)\n');
    fprintf('=========================================================================\n');
    
    addpath(fileparts(mfilename('fullpath')));
    MegaGraph = ApplicationToMegaGraph();
    
    Npop = 100;
    MaxFE = 10000;
    
    % List of Algorithms to Test
    algorithms = {@NSGAII, @NSGAIII, @SPEA2};
    algoNames = {'NSGA-II', 'NSGA-III', 'SPEA2'};
    numAlgos = length(algorithms);
    
    global USE_SPECTRAL_SEEDING
    global CONVERGENCE_TRACE
    
    % Storage for Results
    hv_results = zeros(numAlgos, 2); % Column 1: Random, Column 2: Spectral
    time_results = zeros(numAlgos, 2);
    early_lat_results = zeros(numAlgos, 2); % Latency at 20% of evaluations
    final_lat_results = zeros(numAlgos, 2); % Latency at 100% of evaluations
    
    all_fronts = cell(numAlgos, 2); % To calculate a global reference front for IGD
    traces = cell(numAlgos, 2); % Store convergence traces
    
    for i = 1:numAlgos
        algo = algorithms{i};
        algoName = algoNames{i};
        
        fprintf('\n--- Testing Algorithm: %s ---\n', algoName);
        
        % 1. PURE RANDOM
        fprintf('  Running with Pure Random Seeding...\n');
        USE_SPECTRAL_SEEDING = false;
        CONVERGENCE_TRACE = struct('latency', [], 'thermal', []);
        t1 = tic;
        [~, Obj_r, ~] = platemo('N', Npop, 'maxFE', MaxFE, 'problem', @NocMappingProblem, 'algorithm', algo, 'parameter', MegaGraph, 'save', 0);
        time_results(i, 1) = toc(t1);
        Front_r = GetNonDominated(Obj_r);
        all_fronts{i, 1} = Front_r;
        traces{i, 1} = CONVERGENCE_TRACE.latency;
        
        % 2. SPECTRAL SEEDING
        fprintf('  Running with Spectral Seeding...\n');
        USE_SPECTRAL_SEEDING = true;
        CONVERGENCE_TRACE = struct('latency', [], 'thermal', []);
        t2 = tic;
        [~, Obj_s, ~] = platemo('N', Npop, 'maxFE', MaxFE, 'problem', @NocMappingProblem, 'algorithm', algo, 'parameter', MegaGraph, 'save', 0);
        time_results(i, 2) = toc(t2);
        Front_s = GetNonDominated(Obj_s);
        all_fronts{i, 2} = Front_s;
        traces{i, 2} = CONVERGENCE_TRACE.latency;
    end
    
    % --- POST-PROCESSING (HV and IGD) ---
    fprintf('\nCalculating Metrics (Hypervolume & IGD)...\n');
    
    % Pool all fronts to create a Global Reference Front
    pooled_fronts = [];
    for i = 1:numAlgos
        pooled_fronts = [pooled_fronts; all_fronts{i, 1}; all_fronts{i, 2}];
    end
    RefFront = GetNonDominated(pooled_fronts);
    refPoint = max(RefFront) .* 1.1;
    
    igd_results = zeros(numAlgos, 2);
    
    for i = 1:numAlgos
        % Hypervolume
        hv_results(i, 1) = Calculate2DHV(all_fronts{i, 1}, refPoint);
        hv_results(i, 2) = Calculate2DHV(all_fronts{i, 2}, refPoint);
        
        % IGD
        igd_results(i, 1) = CalculateIGD(all_fronts{i, 1}, RefFront);
        igd_results(i, 2) = CalculateIGD(all_fronts{i, 2}, RefFront);
        
        % Early vs Final Latency (measure of convergence speed)
        early_idx1 = max(1, round(length(traces{i,1}) * 0.2));
        early_idx2 = max(1, round(length(traces{i,2}) * 0.2));
        early_lat_results(i, 1) = traces{i, 1}(early_idx1);
        early_lat_results(i, 2) = traces{i, 2}(early_idx2);
        
        final_lat_results(i, 1) = traces{i, 1}(end);
        final_lat_results(i, 2) = traces{i, 2}(end);
    end
    
    % --- PRINT RESULTS TABLE ---
    fprintf('\n========================================================================================================================================\n');
    fprintf('                                                   ALGORITHMIC COMPARISON RESULTS\n');
    fprintf('========================================================================================================================================\n');
    fprintf('%-12s | %-20s | %-12s | %-12s | %-16s | %-16s | %-10s\n', 'Algorithm', 'Initialization', 'Hypervolume', 'IGD', 'Lat @ 20% FE', 'Lat @ 100% FE', 'Time (s)');
    fprintf('----------------------------------------------------------------------------------------------------------------------------------------\n');
    
    for i = 1:numAlgos
        fprintf('%-12s | %-20s | %-12.4e | %-12.4e | %-16.4f | %-16.4f | %-10.2f\n', algoNames{i}, 'Pure Random', hv_results(i, 1), igd_results(i, 1), early_lat_results(i, 1), final_lat_results(i, 1), time_results(i, 1));
        fprintf('%-12s | %-20s | %-12.4e | %-12.4e | %-16.4f | %-16.4f | %-10.2f\n', '', 'Spectral Seeding', hv_results(i, 2), igd_results(i, 2), early_lat_results(i, 2), final_lat_results(i, 2), time_results(i, 2));
        fprintf('----------------------------------------------------------------------------------------------------------------------------------------\n');
    end
    
    fprintf('Note: Higher Hypervolume is better. Lower IGD is better. Latency @ 20%% shows early convergence speed.\n');
    
    % --- GENERATE CONVERGENCE PLOT ---
    fprintf('\nGenerating Convergence Comparison Plot...\n');
    figure('Name', 'Algorithmic Convergence Comparison', 'Position', [100, 100, 1200, 400]);
    for i = 1:numAlgos
        subplot(1, numAlgos, i);
        plot(traces{i, 1}, 'r-', 'LineWidth', 1.5); hold on;
        plot(traces{i, 2}, 'b-', 'LineWidth', 1.5);
        title(algoNames{i}, 'FontSize', 14);
        xlabel('Evaluation Batches', 'FontSize', 12);
        ylabel('Normalized Latency', 'FontSize', 12);
        grid on;
        legend({'Pure Random', 'Spectral Seeding'}, 'FontSize', 10);
    end
    
    OutDir = fullfile(pwd, 'results', 'figures');
    if ~exist(OutDir, 'dir'); mkdir(OutDir); end
    plotFile = fullfile(OutDir, 'algo_convergence_comparison.png');
    saveas(gcf, plotFile);
    fprintf('Saved Convergence Plot to: %s\n', plotFile);
end

% --- Helper Functions ---
function NonDom = GetNonDominated(Obj)
    if isempty(Obj); NonDom = []; return; end
    N = size(Obj, 1);
    isDominated = false(N, 1);
    for i = 1:N
        for j = 1:N
            if i ~= j
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
    if isempty(Front); hv = 0; return; end
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
    if isempty(Front) || isempty(RefFront)
        igd = Inf;
        return;
    end
    distances = pdist2(RefFront, Front);
    minDistances = min(distances, [], 2);
    igd = mean(minDistances);
end
