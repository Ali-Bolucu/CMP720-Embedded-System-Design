%% VisualizeResults.m
% Generates two figures and saves them to the results folder.
%
%   Figure 1 — Pareto Front scatter plot
%              All solutions in grey, best solution highlighted in red.
%
%   Figure 2 — NoC Grid mapping (9x9)
%              Each cell coloured by application, task ID printed inside.
%
% Based on AA/Visualizer.m reference implementation.
%
% Input:
%   Obj      - [Npop × 2]    Pareto-front objectives  [Latency, FT]
%   BestDec  - [1 × nTasks]  best mapping (task → router ID)
%   BestObj  - [1 × 2]       best individual's objectives
%   MegaGraph - cell array from ApplicationToMegaGraph()
%   OutDir   - output folder path
%   RunLabel - string tag for file names

function VisualizeResults(Obj, BestDec, BestObj, MegaGraph, OutDir, RunLabel)

    if nargin < 5 || isempty(OutDir)
        OutDir = fullfile('results', 'figures');
    end
    if nargin < 6 || isempty(RunLabel)
        RunLabel = datestr(now, 'yyyymmdd_HHMMSS');
    end
    if ~exist(OutDir, 'dir')
        mkdir(OutDir);
    end

    %% ── Parse MegaGraph ─────────────────────────────────────────────────
    nRows        = MegaGraph{2};        % 9
    nCols        = MegaGraph{3};        % 9
    GraphLengths = MegaGraph{8};        % [16 28 53 73]
    nApps        = length(GraphLengths);
    nRoutersFull = nRows * nCols;       % 81

    appEnd   = GraphLengths;
    appStart = [1, GraphLengths(1:end-1) + 1];

    %% ── Sanitise BestDec ─────────────────────────────────────────────────
    % PlatEMO's SBX stores float decision variables internally.
    % Round + clamp to valid router range before building the grid.
    BestDec = round(max(1, min(nRoutersFull, BestDec)));

    %% ── Compute Stats for Best Solution ──────────────────────────────────
    S = MegaGraph{4};
    T = MegaGraph{5};
    W = MegaGraph{6};
    
    [LN, CL] = ind2sub([nRows, nCols], 1:nRoutersFull);
    Pos_Tab = [LN' CL'];
    Dist_Tab = pdist2(Pos_Tab, Pos_Tab, 'cityblock');
    
    % Latency Stats
    sProc = BestDec(S);
    tProc = BestDec(T);
    edgeCost = zeros(1, length(S));
    for edge = 1:length(S)
        edgeCost(edge) = Dist_Tab(sProc(edge), tProc(edge)) * W(edge);
    end
    max_lat = max(edgeCost);
    mean_lat = mean(edgeCost);
    
    % Thermal Hotspot Stats
    taskHeat = zeros(1, length(BestDec));
    for k = 1:length(W)
        taskHeat(S(k)) = taskHeat(S(k)) + W(k);
        taskHeat(T(k)) = taskHeat(T(k)) + W(k);
    end
    alpha = 0.3;
    AdjMatrix = double(Dist_Tab == 1);
    coreHeat = zeros(1, nRoutersFull);
    for t = 1:length(BestDec)
        c = BestDec(t);
        if c >= 1 && c <= nRoutersFull
            coreHeat(c) = coreHeat(c) + taskHeat(t);
        end
    end
    neighborHeat = coreHeat * AdjMatrix;
    coreTemp = coreHeat + alpha * neighborHeat;
    max_temp = max(coreTemp);
    mean_temp = mean(coreTemp);
    
    fprintf('\n[VisualizeResults] ── BEST SOLUTION STATS ──\n');
    fprintf('  Latency : %.2f  (Max: %.0f, Avg: %.2f)\n', BestObj(1), max_lat, mean_lat);
    fprintf('  Thermal : %.2f  (Max: %.2f, Avg: %.2f)\n', BestObj(2), max_temp, mean_temp);
    fprintf('────────────────────────────────────────────\n\n');

    %% ═══════════════════════════════════════════════════════════════════
    %% FIGURE 1 — Pareto Front
    %% ═══════════════════════════════════════════════════════════════════
    fig1 = figure('Name', 'Pareto Front', 'Visible', 'off');

    scatter(Obj(:,1), Obj(:,2), 50, 'k', 'filled');
    hold on;
    scatter(BestObj(1), BestObj(2), 120, 'r', 'p', 'filled');
    hold off;

    xlabel('Latency', 'FontSize', 14);
    ylabel('Thermal Hotspot (TH)', 'FontSize', 14);
    title(sprintf('Pareto Front (NSGA-II) — Run: %s', RunLabel), 'FontSize', 13);
    legend({'Pareto solutions', ...
            sprintf('Best (L=%.2f, TH=%.2f)', BestObj(1), BestObj(2))}, ...
           'Location', 'northeast', 'FontSize', 11);
    grid on;

    pFile = fullfile(OutDir, ['pareto_' RunLabel '.png']);
    exportgraphics(fig1, pFile, 'Resolution', 300, 'BackgroundColor', 'white');
    close(fig1);
    fprintf('[VisualizeResults] Saved Pareto front → %s\n', pFile);

    %% ═══════════════════════════════════════════════════════════════════
    %% FIGURE 2 — NoC Grid Mapping
    %% ═══════════════════════════════════════════════════════════════════
    fig2 = figure('Name', 'NoC Mapping', 'Visible', 'off');

    % Colour palette: one colour per application + white for empty
    Cores = [
        0    0.5  0.5;    % App 1 VOPD   – teal
        0    0.35 0.85;   % App 2 MPEG4  – blue
        0.33 0.65 0.65;   % App 3 VCE    – cyan-green
        0.45 0.55 0.65;   % App 4 WIFIRX – slate
    ];
    emptyColor = [1 1 1];

    % Build PosiNoc: her router için hangi task(lar) atandığını bul
    % PosiNoc{router_idx} = [task_id1, task_id2, ...]  (genellikle tek eleman)
    PosiNoc  = cell(nRoutersFull, 1);
    for taskID = 1:length(BestDec)
        rID = BestDec(taskID);          % already rounded+clamped above
        PosiNoc{rID}(end+1) = taskID;
    end

    % Diagnostic
    nUsed = sum(~cellfun(@isempty, PosiNoc));
    fprintf('[VisualizeResults] Unique routers used: %d / %d   (tasks: %d)\n', ...
            nUsed, nRoutersFull, length(BestDec));

    % MOstraNoC / AA Visualizer ile aynı şekilde yeniden şekillendir.
    % reshape: router ID'ler column-major sıraya göre yerleşir → transpose ile
    % satır/sütun indeksleri düzelir (baseline MOstraNoC permute([2,1,3]) ile aynı).
    numerosCell = reshape(PosiNoc, nRows, nCols);
    numerosCell = numerosCell';

    % App task ranges table
    tarefa_app = zeros(2, nApps);
    for j = 1:nApps
        tarefa_app(1,j) = appStart(j);
        tarefa_app(2,j) = appEnd(j);
    end

    % Cell size
    w = 2;
    h = 2;

    hold on;
    for i = 1:nRows
        for j = 1:nCols
            x = w * (j - 1);
            y = h * (nRows - i);

            tasks = numerosCell{i, j};   % [] = boş, [t1 t2 ...] = dolu

            % Task'ın uygulamasına göre renk belirle (ilk task'a göre)
            Corquadrado = emptyColor;
            if ~isempty(tasks)
                for l = 1:nApps
                    if tasks(1) >= tarefa_app(1,l) && tasks(1) <= tarefa_app(2,l)
                        Corquadrado = Cores(mod(l-1, size(Cores,1))+1, :);
                        break;
                    end
                end
            end

            rectangle('Position', [x y w h], ...
                      'FaceColor', Corquadrado, 'LineWidth', 2);

            if ~isempty(tasks)
                % Tüm task ID'lerini alt alta yaz
                label_tasks = strjoin(arrayfun(@num2str, tasks, 'UniformOutput', false), '/');
                label_weights = strjoin(arrayfun(@(t) sprintf('W:%d', round(taskHeat(t))), tasks, 'UniformOutput', false), '/');
                
                label = {label_tasks, sprintf('(%s)', label_weights)};
                
                text(x + w/2, y + h/2, label, ...
                     'HorizontalAlignment', 'center', ...
                     'FontSize', 7, 'FontWeight', 'bold', 'Color', 'k');
            end
        end
    end
    hold off;

    axis([-1  nCols*w+1  -1  nRows*h+1]);
    axis off;
    title(sprintf('Best Mapping | L=%.2f (Max: %.0f, Avg: %.1f) | TH=%.2f (Max: %.1f, Avg: %.1f)', ...
                  BestObj(1), max_lat, mean_lat, BestObj(2), max_temp, mean_temp), 'FontSize', 11, 'Color', 'k');

    mFile = fullfile(OutDir, ['noc_mapping_' RunLabel '.png']);
    exportgraphics(fig2, mFile, 'Resolution', 300, 'BackgroundColor', 'white');
    close(fig2);
    fprintf('[VisualizeResults] Saved NoC mapping  → %s\n', mFile);

end
