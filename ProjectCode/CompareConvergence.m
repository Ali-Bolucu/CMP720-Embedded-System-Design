%% CompareConvergence.m
% Runs NSGA-II with and without Spectral Seeding to compare convergence speed.

addpath(fileparts(mfilename('fullpath')));

% Configuration
Npop   = 100;
MaxFE  = 10000;
OutDir = fullfile(fileparts(mfilename('fullpath')), 'results', 'figures');
if ~exist(OutDir, 'dir'); mkdir(OutDir); end

RunLabel = datestr(now, 'yyyymmdd_HHMMSS');
logFile = fullfile(OutDir, ['execution_log_' RunLabel '.txt']);
if exist(logFile, 'file'); delete(logFile); end
diary(logFile);
fprintf('--- CompareConvergence Execution Started ---\n');

MegaGraph = ApplicationToMegaGraph();

%% 0.A Spectral Eigenspace Visualization (Diagnostic)
fprintf('\n--- Generating Spectral Clustering Visualization ---\n');
S_edge = MegaGraph{4};
T_edge = MegaGraph{5};
W_edge = MegaGraph{6};
num_tasks = MegaGraph{1};
k_clusters = 4; % We have 4 apps, so k=4 is appropriate

AffinityMatrix = zeros(num_tasks, num_tasks);
for i = 1:length(S_edge)
    src = S_edge(i); tgt = T_edge(i); weight = W_edge(i);
    AffinityMatrix(src, tgt) = AffinityMatrix(src, tgt) + weight;
    AffinityMatrix(tgt, src) = AffinityMatrix(tgt, src) + weight;
end
for i = 1:num_tasks
    AffinityMatrix(i, i) = sum(AffinityMatrix(i, :));
end

D = diag(sum(AffinityMatrix, 2));
D_inv_sqrt = D;
for i = 1:num_tasks
    if D(i,i) > 0; D_inv_sqrt(i,i) = 1 / sqrt(D(i,i)); end
end
L_norm = eye(num_tasks) - D_inv_sqrt * AffinityMatrix * D_inv_sqrt;

[eigenvecs, eigenvals] = eig(L_norm);
eigenvals = diag(eigenvals);
[~, idx] = sort(eigenvals);
idx = idx(2:k_clusters+1);
eigenvecs_selected = eigenvecs(:, idx);

[ClusteredTasks, ~] = kmeans(eigenvecs_selected, k_clusters, 'MaxIter', 300, 'Replicates', 10);
silhouette_vals = silhouette(eigenvecs_selected, ClusteredTasks);
mean_silhouette = mean(silhouette_vals);

fprintf('Mean Silhouette Coefficient: %.4f\n', mean_silhouette);

figure('Visible', 'off', 'Name', 'Spectral Eigenspace', 'Position', [100, 100, 800, 600]);
scatter(eigenvecs_selected(:,1), eigenvecs_selected(:,2), 150, ClusteredTasks, 'filled', 'MarkerEdgeColor', 'k');
colormap('jet');
xlabel('Fiedler Vector 1', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Fiedler Vector 2', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Task Distribution in Spectral Eigenspace\nSilhouette: %.2f', mean_silhouette), 'FontSize', 14);
grid on;
colorbar;

specPlotFile = fullfile(OutDir, ['spectral_eigenspace_' RunLabel '.png']);
exportgraphics(gcf, specPlotFile, 'Resolution', 300, 'BackgroundColor', 'white');
close(gcf);
fprintf('Saved Spectral Eigenspace plot to: %s\n', specPlotFile);

global USE_SPECTRAL_FLAG;
global CONVERGENCE_TRACE;
global PURE_SEED_VALUES;
PURE_SEED_VALUES = struct('Spectral', [0 0], 'Thermal', [0 0]);

global SPECTRAL_SEEDS_EVALS;
SPECTRAL_SEEDS_EVALS = [];

%% 0. Generate Base Random Pool (The "Fairness" Guarantee)
fprintf('\nGenerating strictly identical Base Random Pool (100 individuals)...\n');
global BASE_RANDOM_POOL;
NumTasks = MegaGraph{1};
NumRouters = MegaGraph{2} * MegaGraph{3};
BASE_RANDOM_POOL = zeros(Npop, NumTasks);
for i = 1:Npop
    BASE_RANDOM_POOL(i, :) = randperm(NumRouters, NumTasks);
end

%% 1. Run Pure Random (100% Random)
fprintf('\n======================================================\n');
fprintf('RUNNING PHASE 1: PURE RANDOM (Exploration Only)\n');
fprintf('======================================================\n');
USE_SPECTRAL_FLAG = false;
CONVERGENCE_TRACE = struct('latency', [], 'thermal', []);

% Run PlatEMO
tPhase1 = tic;
[Dec_random, Obj_random, Con_random] = platemo('N', Npop, 'maxFE', MaxFE, 'problem', @NocMappingProblem, 'algorithm', @NSGAII, 'parameter', MegaGraph, 'save', 0);
time_random = toc(tPhase1);
trace_random_lat = CONVERGENCE_TRACE.latency;
trace_random_th  = CONVERGENCE_TRACE.thermal;

%% 2. Run Spectral Seeding (20% Spectral, 80% Random)
fprintf('\n======================================================\n');
fprintf('RUNNING PHASE 2: SPECTRAL SEEDING (20%% Exploitation, 80%% Exploration)\n');
fprintf('======================================================\n');
USE_SPECTRAL_FLAG = true;
CONVERGENCE_TRACE = struct('latency', [], 'thermal', []);

% Run PlatEMO
tPhase2 = tic;
[Dec_spectral, Obj_spectral, Con_spectral] = platemo('N', Npop, 'maxFE', MaxFE, 'problem', @NocMappingProblem, 'algorithm', @NSGAII, 'parameter', MegaGraph, 'save', 0);
time_spectral = toc(tPhase2);
trace_spectral_lat = CONVERGENCE_TRACE.latency;
trace_spectral_th  = CONVERGENCE_TRACE.thermal;

%% 3. Plot Results
fprintf('\nGenerating Convergence Plots...\n');

fig = figure('Name', 'Convergence Comparison', 'Position', [100, 100, 1000, 450], 'Visible', 'off');

% Plot Latency
subplot(1, 2, 1);
plot(1:length(trace_random_lat), trace_random_lat, 'LineWidth', 2, 'Color', [0.8500, 0.3250, 0.0980]);
hold on;
plot(1:length(trace_spectral_lat), trace_spectral_lat, 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410]);
hold off;
title('Latency Convergence', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Evaluations', 'FontSize', 11);
ylabel('Best Latency (Lower is Better)', 'FontSize', 11);
legend('Pure Random', 'Spectral Seeding', 'Location', 'northeast');
grid on;

% Plot Thermal Hotspot
subplot(1, 2, 2);
plot(1:length(trace_random_th), trace_random_th, 'LineWidth', 2, 'Color', [0.8500, 0.3250, 0.0980]);
hold on;
plot(1:length(trace_spectral_th), trace_spectral_th, 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410]);
hold off;
title('Thermal Hotspot Convergence', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Evaluations', 'FontSize', 11);
ylabel('Best Thermal Hotspot (Lower is Better)', 'FontSize', 11);
legend('Pure Random', 'Spectral Seeding', 'Location', 'northeast');
grid on;

% Save Plot
pFile = fullfile(OutDir, ['convergence_' RunLabel '.png']);
exportgraphics(fig, pFile, 'Resolution', 300, 'BackgroundColor', 'white');
close(fig);

fprintf('[CompareConvergence] Saved convergence plot to: %s\n', pFile);

%% 4. Visualize Mapping and Pareto Fronts
fprintf('\nGenerating Pareto Front and Mapping Plots...\n');

% Pure Random Best
BestDec_random = Dec_random(1, :);
BestObj_random = Obj_random(1, :);
VisualizeResults(Obj_random, BestDec_random, BestObj_random, MegaGraph, OutDir, ['Random_' RunLabel]);

% Spectral Seeding Best
BestDec_spectral = Dec_spectral(1, :);
BestObj_spectral = Obj_spectral(1, :);
VisualizeResults(Obj_spectral, BestDec_spectral, BestObj_spectral, MegaGraph, OutDir, ['Spectral_' RunLabel]);

fprintf('[CompareConvergence] Done! All plots saved.\n');

%% 5. Print Comparison Table
fprintf('\n\n==========================================================================================================================================================================\n');
fprintf('                                                             PERFORMANCE & EXPLORATION COMPARISON TABLE\n');
fprintf('==========================================================================================================================================================================\n');
fprintf('%-35s | %-12s | %-12s | %-11s | %-14s | %-14s | %-11s | %-11s | %-12s\n', 'Stage / Approach', 'Best Latency', 'Best Thermal', 'Pareto Size', 'Latency Spread', 'Thermal Spread', 'Hypervolume', 'IGD', 'Exec Time(s)');
fprintf('--------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');

% Initial seeds
best_seed_lat = min(SPECTRAL_SEEDS_EVALS(:, 1));
best_seed_th  = min(SPECTRAL_SEEDS_EVALS(:, 2));

% Extract non-dominated fronts
FrontObj_spectral = GetNonDominated(Obj_spectral);
FrontObj_random   = GetNonDominated(Obj_random);

% Create local Reference Front for IGD (Non-dominated solutions from all pools combined)
RefFront = GetNonDominated([FrontObj_spectral; FrontObj_random]);

% Calculate Reference Point for HV (10% worse than the maximums of both fronts)
refPoint = max(RefFront) .* 1.1;

% Spectral Seeding Final
final_spectral_lat = min(FrontObj_spectral(:, 1));
final_spectral_th  = min(FrontObj_spectral(:, 2));
size_spectral = size(FrontObj_spectral, 1);
spread_spectral_lat = max(FrontObj_spectral(:, 1)) - min(FrontObj_spectral(:, 1));
spread_spectral_th  = max(FrontObj_spectral(:, 2)) - min(FrontObj_spectral(:, 2));
hv_spectral = Calculate2DHV(FrontObj_spectral, refPoint);
igd_spectral = CalculateIGD(FrontObj_spectral, RefFront);

% Pure Random Final
final_random_lat = min(FrontObj_random(:, 1));
final_random_th  = min(FrontObj_random(:, 2));
size_random = size(FrontObj_random, 1);
spread_random_lat = max(FrontObj_random(:, 1)) - min(FrontObj_random(:, 1));
spread_random_th  = max(FrontObj_random(:, 2)) - min(FrontObj_random(:, 2));
hv_random = Calculate2DHV(FrontObj_random, refPoint);
igd_random = CalculateIGD(FrontObj_random, RefFront);

% Print Rows
fprintf('%-35s | %-12.2f | %-12.2f | %-11s | %-14s | %-14s | %-11s | %-11s | %-12s\n', 'Initial 20 Spectral Seeds (Best)', best_seed_lat, best_seed_th, '-', '-', '-', '-', '-', '-');
fprintf('%-35s | %-12.2f | %-12.2f | %-11d | %-14.2f | %-14.2f | %-11.2e | %-11.2e | %-12.2f\n', 'After PlatEMO (Spectral Seeding)', final_spectral_lat, final_spectral_th, size_spectral, spread_spectral_lat, spread_spectral_th, hv_spectral, igd_spectral, time_spectral);
fprintf('%-35s | %-12.2f | %-12.2f | %-11d | %-14.2f | %-14.2f | %-11.2e | %-11.2e | %-12.2f\n', 'After PlatEMO (Pure Random)', final_random_lat, final_random_th, size_random, spread_random_lat, spread_random_th, hv_random, igd_random, time_random);
fprintf('==========================================================================================================================================================================\n');

diary off;

%% ==========================================================================
%% HELPER FUNCTIONS FOR METRICS
%% ==========================================================================

function FrontObj = GetNonDominated(Obj)
    % Extracts the non-dominated solutions (Pareto Front) from a population
    N = size(Obj, 1);
    isDominated = false(N, 1);
    for i = 1:N
        for j = 1:N
            if i ~= j
                % Check if j strictly dominates i
                if all(Obj(j,:) <= Obj(i,:)) && any(Obj(j,:) < Obj(i,:))
                    isDominated(i) = true;
                    break;
                end
            end
        end
    end
    FrontObj = Obj(~isDominated, :);
end

function hv = Calculate2DHV(Front, RefPoint)
    % Simple 2D Hypervolume calculation (minimization)
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
    % IGD is the average minimum distance from points in the True/Reference Front
    % to the nearest point in the obtained Front.
    % Smaller is better.
    
    % Ensure data is present
    if isempty(Front) || isempty(RefFront)
        igd = Inf;
        return;
    end
    
    % pdist2 returns distances [size(RefFront,1) x size(Front,1)]
    % Note: pdist2 might require Statistics and Machine Learning Toolbox
    distances = pdist2(RefFront, Front);
    minDistances = min(distances, [], 2);
    igd = mean(minDistances);
end
