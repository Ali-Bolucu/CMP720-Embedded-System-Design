%% ApplicationGraphToUniGraph.m
% Dört uygulamanın görev grafiklerini tanımlamak ve birleştirmek için kullanılır.

% VOPD (16)        MPEG4 (12)       VCE (25)         WIFIRX (20)
% ├─ Görev 1-16    ├─ Görev 17-28   ├─ Görev 29-53   ├─ Görev 54-73
% └─ Kenarları     └─ Kenarları     └─ Kenarları     └─ Kenarları
%                                 ↓
%                     ┌────────────────────────┐
%                     │  MEGAGraf (73 görev)   │
%                     │  ─ Tüm kenarlar birden │
%                     │  ─ Tüm ağırlıklar      │
%                     └────────────────────────┘

function [MegaGraph] = ApplicationGraphToMegaGraph()
    % Output: [ TasskNumber, SourceTask, TargetTask, EdgeWeight ]

    %% ADIM 1: Uygulama Grafikleri Tanımla
    % APP 1: Video Object Plane Decoder (VOPD)
    apps(1).name = 'VOPD';
    apps(1).nTasks = 16;
    apps(1).source = [1 2 3 4 4 5 6 7 8 8 9 10 11 12 12 12 13 14 15 15 16];
    apps(1).target = [2 3 4 5 16 6 7 8 9 10 10 9 12 6 9 13 14 15 11 13 5];
    apps(1).weight = [70 362 362 362 49 357 353 300 313 500 313 94 16 16 16 16 157 16 16 16 27];
    fprintf('APP 1 - Name: %s, Task number: %d, Edge number: %d\n', apps(1).name, apps(1).nTasks, length(apps(1).source));

    % APP 2: MPEG-4 Video Processing
    apps(2).name = 'MPEG4';
    apps(2).nTasks = 12;
    apps(2).source = [1 1 1 1 1 1 1 2 3 4 5 5 6 6 7 7 8 8 9 9 9 9 10 10 11 12];
    apps(2).target = [2 3 4 5 7 8 10 1 1 1 1 6 5 7 1 6 1 9 8 10 11 12 1 9 9 9];
    apps(2).weight = [64 3 1 20 200 304 11 64 3 1 20 14 14 40 200 40 304 224 224 58 84 167 11 58 84 167];
    fprintf('APP 2 - Name: %s, Task number: %d, Edge number: %d\n', apps(2).name, apps(2).nTasks, length(apps(2).source));
    

    % APP 3: Video Coding Engine (VCE)
    apps(3).name = 'VCE';
    apps(3).nTasks = 25;
    apps(3).source = [1 2 2 3 4 5 6 7 8 8 8 8 9 10 10 11 12 12 13 14 15 15 16 17 18 19 20 22 23 24 25];
    apps(3).target = [2 3 4 4 5 6 18 8 9 12 10 11 12 13 24 10 10 16 14 15 16 17 18 22 19 20 21 23 24 25 9];
    apps(3).weight = [90 90 90 90 30 20 20 8400 2800 2800 2800 5600 2000 4200 4200 1400 30 30 4200 2100 660 660 600 240 620 640 640 240 2210 2280 2280];
    fprintf('APP 3 - Name: %s, Task number: %d, Edge number: %d\n', apps(3).name, apps(3).nTasks, length(apps(3).source));

    % APP 4: WiFi Receiver (WIFIRX)
    apps(4).name = 'WIFIRX';
    apps(4).nTasks = 20;
    apps(4).source = [1 1 2 3 4 4 5 6 7 8 9 9 10 11 12 13 14 15 16 17 17 18 18 18 18 18 18 19 19 19 19 19 19];
    apps(4).target = [2 6 3 4 5 1 6 7 8 9 10 11 11 12 13 14 15 16 17 18 19 12 13 14 15 16 19 1 9 11 17 18 20];
    apps(4).weight = [640 640 640 640 640 1 640 640 512 512 384 384 384 384 72 72 72 108 54 6 54 1 1 1 1 1 4 1 1 1 1 1 54];
    fprintf('APP 4 - Name: %s, Task number: %d, Edge number: %d\n', apps(4).name, apps(4).nTasks, length(apps(4).source));
    

    %% ADIM 2: Grafları Birleştir
    n = length(apps);
    GraphLengths = zeros(1, n); % Her uygulamanın tasklarının başlangıç indeksini tutmak için
    
    % i == 1 durumu: İlklendirmek için
    Task = apps(1).nTasks;
    S = apps(1).source;
    T = apps(1).target;
    W = apps(1).weight;
    GraphLengths(1) = Task; 
    
    % i > 1 durumu: Diğer grafları kaydırma (offset) ile ekle
    for i = 2:n
        S1 = apps(i).source;
        T1 = apps(i).target;
        W1 = apps(i).weight;
        GraphLengths(i) = Task; 
        
        % Düğüm numaralarına o ana kadarki toplam görev sayısını (Task) ekle
        S = cat(2, S, S1 + Task);
        T = cat(2, T, T1 + Task);
        W = cat(2, W, W1);
        
        % Toplam MegaGraf görev sayısını güncelle
        Task = Task + apps(i).nTasks;
    end

    %% ADIM 3: PlatEMO Formatında Final Çıktısı (Cell Array)
    Line = 9;       % NoC ağı 9x9 olacak
    Column = 9;     % NoC ağı 9x9 olacak
    Type = 9;       % OperatorGA dosyasında tanımlı

    CompositeGraph = {Task, Line, Column, S, T, W, Type};
    MegaGraph = [CompositeGraph, {GraphLengths}]; % TODO burayı değiştirdim

    fprintf('MegaGraph is created :\n');
    fprintf(' Toplam Görev: %d ', Task);
    fprintf(' Toplam Kenar: %d\n', length(S));

end
    
