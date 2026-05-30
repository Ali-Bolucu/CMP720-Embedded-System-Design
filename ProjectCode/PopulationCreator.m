function PopulationOut = PopulationCreator(PopSize, NumTasks, NumRouters, AppSizes, useSpectral, S, T, W)
%% PopulationCreator - Initial Population Generator for NoC Task Mapping
% 
% Generates initial population for NoC task-to-router mapping problem.
% Can use either 100% Random Mapping or Spectral Seeding Mapping.
%
% INPUTS:
%   PopSize     - Population size (number of individuals)
%   NumTasks    - Total number of tasks (e.g., 73)
%   NumRouters  - Number of routers = nRows × nCols (e.g., 81 for 9×9)
%   AppSizes    - Vector of task counts per application (e.g., [16 12 25 20])
%   useSpectral - Boolean flag to enable Spectral Seeding
%   S, T, W     - Source, Target, and Weights of communication edges
%--------------------------------------------------------------------------

    % Generate population
    PopulationOut = zeros(PopSize, NumTasks);
    
    global BASE_RANDOM_POOL;
    if isempty(BASE_RANDOM_POOL) || size(BASE_RANDOM_POOL, 1) ~= PopSize || size(BASE_RANDOM_POOL, 2) ~= NumTasks
        disp('[PopulationCreator] BASE_RANDOM_POOL not found or invalid, creating temporary one.');
        BASE_RANDOM_POOL_LOCAL = zeros(PopSize, NumTasks);
        for i = 1:PopSize
            BASE_RANDOM_POOL_LOCAL(i, :) = randperm(NumRouters, NumTasks);
        end
    else
        BASE_RANDOM_POOL_LOCAL = BASE_RANDOM_POOL;
    end
    
    if useSpectral == false
        % PURE RANDOM MAPPING (Maximum Diversity)
        PopulationOut = BASE_RANDOM_POOL_LOCAL;
    else
        % SPECTRAL SEEDING MAPPING
        disp('[PopulationCreator] Running Spectral Clustering for Seeding...');
        
        % 1. Create Affinity Matrix A
        A = zeros(NumTasks, NumTasks);
        for i = 1:length(W)
            A(S(i), T(i)) = A(S(i), T(i)) + W(i);
            A(T(i), S(i)) = A(T(i), S(i)) + W(i); % Symmetric
        end
        
        % 2. Degree Matrix D
        deg = sum(A, 2);
        D = diag(deg);
        
        % 3. Laplacian L (L = D - A as per paper)
        L = D - A;
        
        % 4. Eigendecomposition (V, \Lambda)
        [V, Lambda] = eig(L);
        eigenvalues = diag(Lambda);
        [~, sortIdx] = sort(eigenvalues);
        V_sorted = V(:, sortIdx);
        
        % 5. Normalize rows of V to unit norm
        % Evaluate clusters up to Kmax (e.g. 10)
        Kmax = min(10, NumTasks-1);
        V_k = V_sorted(:, 1:Kmax);
        row_norms = vecnorm(V_k, 2, 2);
        row_norms(row_norms == 0) = 1; % Avoid division by zero
        V_norm = V_k ./ row_norms;
        
        % 6. K-Means & Silhouette to find optimal k*
        best_k = 2;
        best_score = -Inf;
        best_clusters = [];
        opts = statset('Display','off');
        
        for k = 2:Kmax
            % Cluster using the first k eigenvectors
            clusters = kmeans(V_norm(:, 1:k), k, 'Options', opts, 'Replicates', 5);
            % Compute Silhouette score
            sil_scores = silhouette(V_norm(:, 1:k), clusters);
            mean_sil = mean(sil_scores);
            
            if mean_sil > best_score
                best_score = mean_sil;
                best_k = k;
                best_clusters = clusters;
            end
        end
        
        fprintf('[PopulationCreator] Optimal clusters k* = %d (Silhouette: %.3f)\n', best_k, best_score);
        
        % 7. Map clusters to NoC (Seeding)
        nRows = ceil(sqrt(NumRouters));
        nCols = floor(sqrt(NumRouters));
        [LN, CL] = ind2sub([nRows, nCols], 1:NumRouters);
        PosMesh = [LN', CL'];
        DistMatrix = pdist2(PosMesh, PosMesh, 'cityblock');
        
        for individual = 1:PopSize
            if individual <= round(0.20 * PopSize)
                % --- 20% SPECTRAL SEEDING ---
                meshGrid = zeros(nRows, nCols);
                clusterOrder = randperm(best_k);
                
                for cIdx = 1:best_k
                    cID = clusterOrder(cIdx);
                    tasksInCluster = find(best_clusters == cID)';
                    tasksInCluster = tasksInCluster(randperm(length(tasksInCluster))); % Shuffle
                    
                    for taskPos = 1:length(tasksInCluster)
                        taskID = tasksInCluster(taskPos);
                        
                        if taskPos == 1
                            % Place first task of cluster randomly
                            placed = false;
                            while ~placed
                                row = randi(nRows);
                                col = randi(nCols);
                                if meshGrid(row, col) == 0
                                    meshGrid(row, col) = taskID;
                                    lastRow = row;
                                    lastCol = col;
                                    placed = true;
                                end
                            end
                        else
                            % Fuzzy spatial locality for the rest of the cluster
                            [newRow, newCol] = FindNearestFreeRouter(meshGrid, lastRow, lastCol, DistMatrix, nRows, nCols);
                            meshGrid(newRow, newCol) = taskID;
                            lastRow = newRow;
                            lastCol = newCol;
                        end
                    end
                end
                
                % Convert 2D grid to 1D mapping vector
                mapping1D = meshGrid(:)';
                indices = find(mapping1D);
                individual_mapping = zeros(1, NumTasks);
                individual_mapping(mapping1D(indices)) = indices;
                PopulationOut(individual, :) = individual_mapping;
                
                % Evaluate this Spectral seed without affecting convergence trace
                global CONVERGENCE_TRACE PURE_SEED_VALUES SPECTRAL_SEEDS_EVALS;
                temp_trace = CONVERGENCE_TRACE;
                CONVERGENCE_TRACE = []; % Disable trace
                objVal = NocObjectiveEvaluator(individual_mapping, DistMatrix, NumRouters, S, T, W, nRows, nCols);
                CONVERGENCE_TRACE = temp_trace; % Restore trace
                
                SPECTRAL_SEEDS_EVALS(individual, :) = objVal;
                
                % Save the BEST (minimum) Latency and Thermal among the 20 spectral seeds
                if individual == 1
                    PURE_SEED_VALUES.Spectral = objVal;
                else
                    if objVal(1) < PURE_SEED_VALUES.Spectral(1)
                        PURE_SEED_VALUES.Spectral(1) = objVal(1);
                    end
                    if objVal(2) < PURE_SEED_VALUES.Spectral(2)
                        PURE_SEED_VALUES.Spectral(2) = objVal(2);
                    end
                end
                
            else
                % --- 80% PURE RANDOM (Exploration) ---
                PopulationOut(individual, :) = BASE_RANDOM_POOL_LOCAL(individual, :);
            end
        end
    end
end

%% ════════════════════════════════════════════════════════════════════════
%% HELPER FUNCTION: Find Nearest Free Router
%% ════════════════════════════════════════════════════════════════════════

function [newRow, newCol] = FindNearestFreeRouter(meshGrid, currentRow, currentCol, ...
                                                   distMatrix, nRows, nCols)
    % Convert current position to 1D index
    current1D = sub2ind([nRows, nCols], currentRow, currentCol);
    distances = distMatrix(current1D, :);
    
    % Find free routers (where meshGrid == 0)
    meshGrid1D = meshGrid(:)';
    freeRouters = find(meshGrid1D == 0);
    
    % Sort free routers by distance
    [~, sortIdx] = sort(distances(freeRouters));
    
    numFree = length(freeRouters);
    % Fuzzy spatial locality: 10% chance to jump randomly, 90% chance to pick from top 4
    if rand() < 0.10
        chosenIdx = randi(numFree);
    else
        K = min(4, numFree);
        topK_indices = sortIdx(1:K);
        chosenIdx = topK_indices(randi(K));
    end
    
    nearest1D = freeRouters(chosenIdx);
    [newRow, newCol] = ind2sub([nRows, nCols], nearest1D);
end
