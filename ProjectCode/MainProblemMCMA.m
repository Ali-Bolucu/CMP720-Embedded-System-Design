%% MainProblemMCMA.m
% MCMA (Many-Core Mapping Optimization)

% Ensure current folder and subfolders are on the MATLAB path
currentDir = fileparts(mfilename('fullpath'));
addpath(currentDir);

%% ADIM 1: Uygulama graflarından PletEmo için tek bir megagraf oluştur
fprintf('STEP1 ==============================\n');
MegaGraph= ApplicationGraphToMegaGraph();
fprintf('====================================\n');



%% ADIM 2: PletEmo'yu çalıştır
fprintf('STEP2 ==============================\n');

for t = 1 : 1
    disp(t);
    %[Dec,Obj,Con] = platemo('problem',@ManyCoreMAV1,'algorithm',@NSGAII,'parameter',{nTask,Line,Column,S,T,P},'save', 1);
    [Dec,Obj,Con] = platemo('N',100,'problem',@MCMAOneTask,'algorithm',@NSGAII,'parameter',MegaGraph,'save', 1);

    % Optimizasyon sonuçları
    AllDecisions = Dec;     % Tüm çözümlerin görev atamaları
    AllObjectives = Obj;    % Tüm çözümlerin amaç değerleri (Latency, FT)
    AllConstraints = Con;   % Kısıtlamalar (genellikle boş)
end

%cd(currentDir);
fprintf('====================================\n');


%% ADIM 3: Sonuçları Görselleştir
fprintf('====================================\n');

% 1. Çalışma dizininde 'Output' klasörü oluştur
outputDir = fullfile(pwd, 'Output');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% 2. Pareto Front'u çiz ve kaydet
Visualizer.plotPareto(AllObjectives, outputDir);

% 3. Çizim için en iyi çözümü belirle
% Pareto cephesindeki çözümlerden Gecikmesi (Latency) en düşük olanı seçiyoruz.
% (Alternatif olarak orijine en yakın olan da seçilebilir)
[~, bestIdx] = min(AllObjectives(:, 1)); 
bestDecision = AllDecisions(bestIdx, :);

% 4. NoC haritasını çiz ve kaydet
Visualizer.plotNoC(bestDecision, MegaGraph, outputDir);

fprintf('====================================\n');

