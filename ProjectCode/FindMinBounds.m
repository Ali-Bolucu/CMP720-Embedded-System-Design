%% FindMinBounds.m
% Automatically finds independent ideal (minimum) bounds for each metric
% (max_latency, mean_latency, max_temp, mean_temp) using Simulated Annealing.
% Used for identifying the 0-point for Min-Max normalization.

function FindMinBounds()
    fprintf('=== Normalization Bounds Finder (Simulated Annealing - MINIMIZATION) ===\n');
    
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
    
    % 1. Minimize Max Latency
    fprintf('\nMinimizing Max Latency...\n');
    best_max_lat = RunSAMin(@(m) EvalObj(m, 1, Dist_Tab, S, T, W, nCores, taskHeat, AdjMatrix), nTasks, nCores, maxIter, initialTemp, coolingRate);
    fprintf('>>> Best Max Latency Found: %.2f\n', best_max_lat);
    
    % 2. Minimize Mean Latency
    fprintf('\nMinimizing Mean Latency...\n');
    best_mean_lat = RunSAMin(@(m) EvalObj(m, 2, Dist_Tab, S, T, W, nCores, taskHeat, AdjMatrix), nTasks, nCores, maxIter, initialTemp, coolingRate);
    fprintf('>>> Best Mean Latency Found: %.2f\n', best_mean_lat);
    
    % 3. Minimize Max Temp
    fprintf('\nMinimizing Max Thermal Hotspot...\n');
    best_max_temp = RunSAMin(@(m) EvalObj(m, 3, Dist_Tab, S, T, W, nCores, taskHeat, AdjMatrix), nTasks, nCores, maxIter, initialTemp, coolingRate);
    fprintf('>>> Best Max Thermal Found: %.2f\n', best_max_temp);
    
    % 4. Minimize Mean Temp
    fprintf('\nMinimizing Mean Thermal (Overall Heating)...\n');
    best_mean_temp = RunSAMin(@(m) EvalObj(m, 4, Dist_Tab, S, T, W, nCores, taskHeat, AdjMatrix), nTasks, nCores, maxIter, initialTemp, coolingRate);
    fprintf('>>> Best Mean Thermal Found: %.2f\n', best_mean_temp);
    
    fprintf('\n========================================================\n');
    fprintf('                   FINAL BOUNDS (IDEAL POINTS)          \n');
    fprintf('========================================================\n');
    fprintf('BOUNDS.min_maxLat  = %.2f;\n', best_max_lat);
    fprintf('BOUNDS.min_meanLat = %.2f;\n', best_mean_lat);
    fprintf('BOUNDS.min_maxTemp = %.2f;\n', best_max_temp);
    fprintf('BOUNDS.min_meanTemp= %.2f;\n', best_mean_temp);
    fprintf('========================================================\n');
end

%% Simulated Annealing Core (Minimization)
function bestCost = RunSAMin(objFunc, nTasks, nCores, maxIter, initialTemp, coolingRate)
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
        
        % Acceptance Probability (Minimization)
        if newCost < currentCost
            currentMapping = newMapping;
            currentCost = newCost;
            if currentCost < bestCost
                bestCost = currentCost;
            end
        else
            prob = exp((currentCost - newCost) / Temp);
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
