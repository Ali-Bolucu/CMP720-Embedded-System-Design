%% GeneratePaperVisuals.m
% Bu script, makalenin "Methodology" bölümünde kullanılmak üzere
% yüksek kaliteli, akademik görseller (Affinity Matrix ve Fiedler Vectors) üretir.
% Çıktılar 'Extended_Project_Report/figures' klasörüne kaydedilir.

clc; clear; close all;

%% 1. Veri Yükleme ve Affinity Matrix Oluşturma
fprintf('Generating MegaGraph...\n');
MegaGraph = ApplicationToMegaGraph();

NumTasks = MegaGraph{1};
S = MegaGraph{4};
T = MegaGraph{5};
W = MegaGraph{6};
AppLengths = MegaGraph{8};

% Simetrik Affinity Matrisi A
A = zeros(NumTasks, NumTasks);
for i = 1:length(W)
    A(S(i), T(i)) = A(S(i), T(i)) + W(i);
    A(T(i), S(i)) = A(T(i), S(i)) + W(i); % Symmetric
end

%% 2. Figure 1: Block-Diagonal Affinity Matrix Heatmap
fprintf('Plotting Affinity Matrix...\n');
fig1 = figure('Name', 'Unified Affinity Matrix', 'Position', [100, 100, 600, 600]);

% Ağırlıkları daha iyi görselleştirmek için log-scale kullan
A_vis = log10(A + 1); 

imagesc(A_vis);
colormap(flipud(hot)); % Hot colormap ters çevrilerek arka plan beyaz, bağlar koyu kırmızı
colorbar('Ticks', [0, max(A_vis(:))], 'TickLabels', {'0', 'Max Bandwidth'});

hold on;
% Uygulama sınırlarını çiz (Block-Diagonal yapı)
app_bounds = [0, AppLengths];
for i = 2:length(app_bounds)-1
    xline(app_bounds(i)+0.5, 'k--', 'LineWidth', 1.5);
    yline(app_bounds(i)+0.5, 'k--', 'LineWidth', 1.5);
end

% Uygulama isimlerini ekle
labels = {'VOPD', 'MPEG4', 'VCE', 'WIFIRX'};
for i = 1:length(labels)
    mid_pos = (app_bounds(i) + app_bounds(i+1)) / 2;
    text(mid_pos, mid_pos, labels{i}, 'HorizontalAlignment', 'center', ...
        'Color', 'blue', 'FontWeight', 'bold', 'FontSize', 12);
end

title('Unified Block-Diagonal Affinity Matrix');
xlabel('Destination Task ID');
ylabel('Source Task ID');
axis square;
set(gca, 'FontSize', 12, 'LineWidth', 1);

%% 3. Laplacian ve Spektral Analiz
fprintf('Computing Laplacian and Spectral Clustering...\n');
deg = sum(A, 2);
D = diag(deg);
L = D - A;

[V, Lambda] = eig(L);
eigenvalues = diag(Lambda);
[~, sortIdx] = sort(eigenvalues);
V_sorted = V(:, sortIdx);

% Kmax'a kadar olan vektörleri al
Kmax = min(10, NumTasks-1);
V_k = V_sorted(:, 1:Kmax);
row_norms = vecnorm(V_k, 2, 2);
row_norms(row_norms == 0) = 1; 
V_norm = V_k ./ row_norms;

% K-Means ile en iyi kümeyi bul
best_k = 2;
best_score = -Inf;
best_clusters = [];
opts = statset('Display','off');

all_scores = zeros(1, Kmax);

for k = 2:Kmax
    clusters = kmeans(V_norm(:, 1:k), k, 'Options', opts, 'Replicates', 5);
    sil_scores = silhouette(V_norm(:, 1:k), clusters);
    mean_sil = mean(sil_scores);
    all_scores(k) = mean_sil;
    if mean_sil > best_score
        best_score = mean_sil;
        best_k = k;
        best_clusters = clusters;
    end
end

fprintf('Optimal clusters found: k = %d\n', best_k);

%% 4. Figure 2: Fiedler Vector Scatter Plot (Spectral Space)
fprintf('Plotting Spectral Space (Fiedler Vectors)...\n');
fig2 = figure('Name', 'Spectral Space Representation', 'Position', [650, 100, 800, 500]);

% Kümelerdeki görev sayılarını hesapla
cluster_counts = groupcounts(best_clusters);

% 2., 3. ve 4. Özvektörleri (Fiedler) kullanarak 3 Boyutlu çizim
scatter3(V_norm(:,2), V_norm(:,3), V_norm(:,4), 100, best_clusters, 'filled', 'MarkerEdgeColor', 'k');
colormap(lines(best_k));
colorbar('Ticks', 1:best_k, 'TickLabels', arrayfun(@(x) sprintf('C%d (n=%d)', x, cluster_counts(x)), 1:best_k, 'UniformOutput', false));

title(sprintf('3D Spectral Clustering in Eigenspace (k*=%d)', best_k));
xlabel('Fiedler Vector 1 ($V_2$)', 'Interpreter', 'latex');
ylabel('Fiedler Vector 2 ($V_3$)', 'Interpreter', 'latex');
zlabel('Fiedler Vector 3 ($V_4$)', 'Interpreter', 'latex');
view(45, 30); % 3D açıyı ayarla
grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1);

%% 4.1 Figure 3: Silhouette Score Evaluation Plot
fprintf('Plotting Silhouette Evaluation...\n');
fig3 = figure('Name', 'Silhouette Score Evaluation', 'Position', [200, 200, 600, 400]);
plot(2:Kmax, all_scores(2:end), '-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
hold on;
plot(best_k, best_score, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r'); % Optimum noktayı yıldızla işaretle
text(best_k, best_score + 0.05, sprintf('Optimal $k^* = %d$', best_k), 'Interpreter', 'latex', 'FontSize', 12, 'HorizontalAlignment', 'center', 'Color', 'red');

title('Dynamic Cluster Determination via Silhouette Coefficient');
xlabel('Number of Clusters ($k$)', 'Interpreter', 'latex');
ylabel('Average Silhouette Score', 'Interpreter', 'latex');
grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1);

%% 5. Dosyaları Kaydet (PDF formatında)
save_dir = 'Extended_Project_Report/figures';
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

fprintf('Saving figures to %s...\n', save_dir);
exportgraphics(fig1, fullfile(save_dir, 'fig_affinity_matrix.pdf'), 'ContentType', 'vector');
exportgraphics(fig1, fullfile(save_dir, 'fig_affinity_matrix.png'), 'Resolution', 300);

exportgraphics(fig2, fullfile(save_dir, 'fig_spectral_space.pdf'), 'ContentType', 'vector');
exportgraphics(fig2, fullfile(save_dir, 'fig_spectral_space.png'), 'Resolution', 300);

exportgraphics(fig3, fullfile(save_dir, 'fig_silhouette_evaluation.pdf'), 'ContentType', 'vector');
exportgraphics(fig3, fullfile(save_dir, 'fig_silhouette_evaluation.png'), 'Resolution', 300);

fprintf('Done! Visuals generated successfully.\n');
