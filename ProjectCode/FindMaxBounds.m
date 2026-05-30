%% FindMaxBounds.m
% Automatically finds independent nadir (maximum) bounds for each metric
% (max_latency, mean_latency, max_temp, mean_temp) using Simulated Annealing.
% Used for identifying the 100-point for Min-Max normalization.

function FindMaxBounds()
    fprintf('=== Normalization Bounds Finder (Simulated Annealing - MAXIMIZATION) ===\n');
    
    % Load Graph Data
    addpath(fileparts(mfilename('fullpath')));
    MegaGraph = ApplicationToMegaGraph();
    nTasks = MegaGraph{1};
    nR = MegaGraph{2};
    nC = MegaGraph{3};
    S = MegaGraph{4};
    T = MegaGraph{5};
    W = MegaGraph{6};
    
    nCores = nR * nC;
    
    % Precompute Distance and Adjacency Matrices
    [LN, CL] = ind2sub([nR, nC], 1:nCores);
    Pos_Tab = [LN' CL'];
    Dist_Tab = pdist2(Pos_Tab, Pos_Tab, 'cityblock');
    AdjMatrix = double(Dist_Tab == 1);
    
    % Precompute Task Heat (Constant)
    taskHeat = zeros(1, nTasks);
    for k = 1:length(W)
        taskHeat(S(k)) = taskHeat(S(k)) + W(k);
        taskHeat(T(k)) = taskHeat(T(k)) + W(k);
    end
    
    % Setup SA Parameters
    maxIter = 50000;
    initialTemp = 5000;
    coolingRate = 0.999;
    
    % 1. Maximize Max Latency
    fprintf('\nMaximizing Max Latency...\n');
    best_max_lat = RunSAMax(@(m) EvalObj(m, 1, Dist_Tab, S, T, W, nCores, taskHeat, AdjMatrix), nTasks, nCores, maxIter, initialTemp, coolingRate);
    fprintf('>>> Best Max Latency Found: %.2f\n', best_max_lat);
    
    % 2. Maximize Mean Latency
    fprintf('\nMaximizing Mean Latency...\n');
    best_mean_lat = RunSAMax(@(m) EvalObj(m, 2, Dist_Tab, S, T, W, nCores, taskHeat, AdjMatrix), nTasks, nCores, maxIter, initialTemp, coolingRate);
    fprintf('>>> Best Mean Latency Found: %.2f\n', best_mean_lat);
    
    % 3. Maximize Max Temp
    fprintf('\nMaximizing Max Thermal Hotspot...\n');
    best_max_temp = RunSAMax(@(m) EvalObj(m, 3, Dist_Tab, S, T, W, nCores, taskHeat, AdjMatrix), nTasks, nCores, maxIter, initialTemp, coolingRate);
    fprintf('>>> Best Max Thermal Found: %.2f\n', best_max_temp);
    
    % 4. Maximize Mean Temp
    fprintf('\nMaximizing Mean Thermal (Overall Heating)...\n');
    best_mean_temp = RunSAMax(@(m) EvalObj(m, 4, Dist_Tab, S, T, W, nCores, taskHeat, AdjMatrix), nTasks, nCores, maxIter, initialTemp, coolingRate);
    fprintf('>>> Best Mean Thermal Found: %.2f\n', best_mean_temp);
    
    fprintf('\n========================================================\n');
    fprintf('                   FINAL BOUNDS (NADIR POINTS)          \n');
    fprintf('========================================================\n');
    fprintf('BOUNDS.max_maxLat  = %.2f;\n', best_max_lat);
    fprintf('BOUNDS.max_meanLat = %.2f;\n', best_mean_lat);
    fprintf('BOUNDS.max_maxTemp = %.2f;\n', best_max_temp);
    fprintf('BOUNDS.max_meanTemp= %.2f;\n', best_mean_temp);
    fprintf('========================================================\n');
end

%% Simulated Annealing Core (Maximization)
function bestCost = RunSAMax(objFunc, nTasks, nCores, maxIter, initialTemp, coolingRate)
    currentMapping = randperm(nCores, nTasks);
    currentCost = objFunc(currentMapping);
    
    bestCost = currentCost;
    Temp = initialTemp;
    
    for i = 1:maxIter
        % Generate Neighbor
        newMapping = currentMapping;
        if rand() < 0.5
            idx = randperm(nTasks, 2);
            temp = newMapping(idx(1));
            newMapping(idx(1)) = newMapping(idx(2));
            newMapping(idx(2)) = temp;
        else
            mappedCores = newMapping;
            emptyCores = setdiff(1:nCores, mappedCores);
            if ~isempty(emptyCores)
                tIdx = randi(nTasks);
                eIdx = randi(length(emptyCores));
                newMapping(tIdx) = emptyCores(eIdx);
            end
        end
        
        newCost = objFunc(newMapping);
        
        % Acceptance Probability (Maximization)
        if newCost > currentCost
            currentMapping = newMapping;
            currentCost = newCost;
            if currentCost > bestCost
                bestCost = currentCost;
            end
        else
            prob = exp((newCost - currentCost) / Temp);
            if rand() < prob
                currentMapping = newMapping;
                currentCost = newCost;
            end
        end
        Temp = Temp * coolingRate;
    end
end

%% Evaluation Helper
function val = EvalObj(mapping, targetMetric, Dist_Tab, S, T, W, nCores, taskHeat, AdjMatrix)
    if targetMetric == 1 || targetMetric == 2
        numEdges = length(S);
        edgeCost = zeros(1, numEdges);
        for edge = 1:numEdges
            src = mapping(S(edge));
            tgt = mapping(T(edge));
            edgeCost(edge) = Dist_Tab(src, tgt) * W(edge);
        end
        if targetMetric == 1; val = max(edgeCost); else; val = mean(edgeCost); end
    elseif targetMetric == 3 || targetMetric == 4
        coreHeat = zeros(1, nCores);
        for t = 1:length(mapping)
            c = mapping(t);
            coreHeat(c) = taskHeat(t);
        end
        alpha = 0.3;
        neighborHeat = coreHeat * AdjMatrix;
        coreTemp = coreHeat + alpha * neighborHeat;
        if targetMetric == 3; val = max(coreTemp); else; val = mean(coreTemp); end
    end
end
