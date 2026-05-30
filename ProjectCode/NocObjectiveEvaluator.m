function PopObj = NoCObjectiveEvaluator(PopDec, Dist_Tab, NumCores, S, T, W, nR, nC)
    %% NoCObjectiveEvaluator - Multi-Objective NoC Task Mapping Fitness Evaluator
    % 
    % Evaluates two objectives for task-to-router assignments on a 2D mesh NoC:
    %   Objective 1: Communication Latency (maximum edge path length × weight)
    %   Objective 2: Thermal Hotspot Cost (heat accumulation across cores)
    %
    % INPUTS:
    %   PopDec     [L×D] - Population decisions: each row = task-to-router assignment
    %   Dist_Tab   [N×N] - Manhattan distance matrix between all cores (N = nR×nC)
    %   NumCores   [scalar] - Total number of cores (nR × nC)
    %   S          [1×E] - Source task indices (edges)
    %   T          [1×E] - Target task indices (edges)
    %   W          [1×E] - Edge weights (communication volume)
    %   nR         [scalar] - NoC mesh rows
    %   nC         [scalar] - NoC mesh columns
    %
    % OUTPUTS:
    %   PopObj     [L×2] - Two objectives for each individual: [Latency, ThermalCost]
    %
    % Created: 2026 (Refactored from MCMACustoV4.m)
    % Original: Manoel Aranda de Almeida 25/05/2023
    %% ========================================================================
    
    [L, D] = size(PopDec);  % L = population size, D = number of tasks
    
    % =========================================================================
    % REPAIR CONTINUOUS / INVALID MAPPINGS (for algorithms like MOEA/D)
    % Some algorithms produce continuous (float) or duplicate values.
    % We must convert them back to unique valid integer core IDs (1 to nR*nC).
    % =========================================================================
    nCores = nR * nC;
    for i = 1:L
        row = PopDec(i, :);
        % If not integers, or has duplicates, or out of bounds
        if any(row ~= round(row)) || length(unique(round(row))) < D || any(row < 1) || any(row > nCores)
            row = round(row);
            row = max(1, min(nCores, row)); % Clamp bounds
            [u, idx] = unique(row, 'stable');
            if length(u) < D
                % Fill duplicates with random unused cores
                available = setdiff(1:nCores, u);
                missing_count = D - length(u);
                if missing_count > length(available)
                    missing_count = length(available); % Safety
                end
                avail_shuffled = available(randperm(length(available)));
                
                new_row = zeros(1, D);
                new_row(idx) = u;
                dup_idx = setdiff(1:D, idx);
                new_row(dup_idx(1:missing_count)) = avail_shuffled(1:missing_count);
                row = new_row;
            end
            PopDec(i, :) = row;
        end
    end
    
    %% Objective 1: Communication Latency
    % Calculate maximum weighted distance for each individual's mapping
    latency = EvaluateCommunicationLatency(PopDec, Dist_Tab, S, T, W);
    
    %% Objective 2: Thermal Hotspot Cost
    % Calculate thermal stress on NoC cores
    thermalCost = EvaluateThermalHotspot(PopDec, nR, nC, Dist_Tab, S, T, W);
    
    
    %% Combine objectives
    PopObj = [latency thermalCost];
    
    %% Convergence Tracker
    global CONVERGENCE_TRACE;
    if isstruct(CONVERGENCE_TRACE)
        % Find best in this evaluation batch
        best_lat = min(latency);
        best_th  = min(thermalCost);
        
        if isempty(CONVERGENCE_TRACE.latency)
            CONVERGENCE_TRACE.latency = best_lat;
            CONVERGENCE_TRACE.thermal = best_th;
        else
            % Track the absolute best seen so far across generations
            CONVERGENCE_TRACE.latency(end+1) = min(CONVERGENCE_TRACE.latency(end), best_lat);
            CONVERGENCE_TRACE.thermal(end+1) = min(CONVERGENCE_TRACE.thermal(end), best_th);
        end
    end
    
end

%% ==========================================================================
%% HELPER FUNCTION 1: Communication Latency Evaluator
%% ==========================================================================

function latency = EvaluateCommunicationLatency(PopDec, Dist_Tab, S, T, W)
    %% EvaluateCommunicationLatency
    % Computes maximum communication latency (bottleneck edge) for each individual
    %
    % For each edge (source→target):
    %   - Get assigned cores: sProc, tProc
    %   - Get path distance from Dist_Tab
    %   - Multiply by edge weight: cost = distance × weight
    % Objective = max cost across all edges (bottleneck)
    
    [L, ~] = size(PopDec);
    numEdges = length(S);
    
    % Extract source and target processors for all edges
    sProc = PopDec(:, S);  % [L × numEdges]
    tProc = PopDec(:, T);  % [L × numEdges]
    
    % Calculate distance for each edge in each individual
    edgeCost = zeros(L, numEdges);
    
    for individual = 1:L
        for edge = 1:numEdges
            src = sProc(individual, edge);
            tgt = tProc(individual, edge);
            edgeCost(individual, edge) = Dist_Tab(src, tgt) * W(edge);
        end
    end
    
    % Latency = Hybrid metric (Weighted sum of Maximum bottleneck and Average cost)
    max_latency = max(edgeCost, [], 2);   % [L × 1]
    mean_latency = mean(edgeCost, 2);     % [L × 1]
    
    % --- NORMALIZATION BOUNDS ---
    % Min bounds found via Simulated Annealing (Ideal Points)
    min_maxLat  = 8400.00;
    min_meanLat = 699.46;
    % Max bounds found via SA (Nadir Points)
    max_maxLat  = 134400.00;
    max_meanLat = 7536.07;
    
    % Min-Max Scaling (0 to 1)
    norm_max_latency = (max_latency - min_maxLat) ./ (max_maxLat - min_maxLat);
    norm_mean_latency = (mean_latency - min_meanLat) ./ (max_meanLat - min_meanLat);
    
    % Clamp values just in case a new mapping goes slightly out of bounds
    norm_max_latency = max(0, min(1, norm_max_latency));
    norm_mean_latency = max(0, min(1, norm_mean_latency));
    
    % Weights for the hybrid objective (e.g., 70% Maximum, 30% Average)
    w_max  = 0.7;
    w_mean = 0.3;
    
    latency = (w_max * norm_max_latency) + (w_mean * norm_mean_latency);
    
end

%% ==========================================================================
%% HELPER FUNCTION 2: Thermal Hotspot Cost Evaluator
%% ==========================================================================

function tCost = EvaluateThermalHotspot(PopDec, nR, nC, Dist_Tab, S, T, W)
    %% EvaluateThermalHotspot
    % Calculates the thermal profile of cores based on task communication
    % volumes and returns the thermal hotspot cost (peak temperature).
    
    nSolutions = size(PopDec, 1);
    nCores     = nR * nC;
    tCost      = zeros(nSolutions, 1);
    
    % Step 1: Calculate total communication volume for each task
    nTasks = size(PopDec, 2);
    taskHeat = zeros(1, nTasks);
    for k = 1:length(W)
        taskHeat(S(k)) = taskHeat(S(k)) + W(k);
        taskHeat(T(k)) = taskHeat(T(k)) + W(k);
    end
    
    % Heat dissipation coefficient (ratio of heat passed to neighbors)
    alpha = 0.3; 
    
    % Adjacency matrix (neighbors at distance == 1)
    % Used for vectorized computation of heat dissipation
    AdjMatrix = double(Dist_Tab == 1); 
    
    for i = 1:nSolutions
        mapping = PopDec(i, :);
        coreHeat = zeros(1, nCores);
        
        % Step 2: Calculate "specific heat" for each core
        for t = 1:length(mapping)
            c = mapping(t);
            if c >= 1 && c <= nCores
                coreHeat(c) = coreHeat(c) + taskHeat(t);
            end
        end
        
        % Step 3: Apply thermal dissipation via matrix multiplication
        % neighborHeat: total heat received from neighboring cores
        neighborHeat = coreHeat * AdjMatrix;
        coreTemp = coreHeat + alpha * neighborHeat;
        
        % Step 4: Calculate hotspot (peak) and overall heating (mean)
        maxTemp  = max(coreTemp);
        meanTemp = mean(coreTemp);
        
        % --- NORMALIZATION BOUNDS ---
        min_maxTemp  = 22400.00;
        min_meanTemp = 2930.75;
        max_maxTemp  = 33836.00;
        max_meanTemp = 3561.73;
        
        % Min-Max Scaling (0 to 1)
        norm_maxTemp = (maxTemp - min_maxTemp) / (max_maxTemp - min_maxTemp);
        norm_meanTemp = (meanTemp - min_meanTemp) / (max_meanTemp - min_meanTemp);
        
        % Clamp values
        norm_maxTemp = max(0, min(1, norm_maxTemp));
        norm_meanTemp = max(0, min(1, norm_meanTemp));
        
        % Objective: Minimize both regional hotspot and overall chip heating
        % Weighted sum: 70% Hotspot, 30% Overall Heating
        weightHotspot = 0.7;
        weightMean    = 0.3;
        
        tCost(i) = weightHotspot * norm_maxTemp + weightMean * norm_meanTemp;
    end
end