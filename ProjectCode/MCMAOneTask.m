%% MCMAOneTask.m
% MCMAOneTask: Many-Core Mapping Optimization for a Single Task Graph

%% AKIŞ Şeması
%   MainProblemMCMA.m
%       ↓
%       "platemo(..., 'problem', @MCMAOneTask, ...)"
%       ↓
%   PlatEMO şuna bakıyor:
%   "Tamam, bu problem sınıfı ne istiyormuş?"
%       ↓
%   MCMAOneTask.m
%       ├─ Setting() → Problem boyutunu tanımla
%       ├─ Initialization() → İlk popülasyonu oluştur
%       └─ CalObj() → Her çözümü değerlendir
%       ↓
%   NSGA-II Algoritması çalışıyor:
%       ├─ Initialization() → 100 rastgele çözüm
%       ├─ CalObj() → Her çözümün iyiliğini ölç
%       ├─ Mutasyon / Seçim
%       ├─ CalObj() → Yeniden değerlendir
%       ├─ Tekrar... (200 nesil)
%       └─ Best çözümleri döndür
%       ↓
%   [Dec, Obj, Con]



classdef MCMAOneTask < PROBLEM % Inherit from PROBLEM class
    
    % MegaGraf objesinin parametreleri
    properties(Access = private)
        nTask   % Toplam görev sayısı 
        Line    % NoC ağı row sayısı
        Column  % NoC ağı column sayısı
        S       % Source görevler (MegaGraf'taki kaynak düğümler)
        T       % Target görevler (MegaGraf'taki hedef düğümler)
        W       % Edge ağırlıkları (MegaGraf'taki görevler arasındaki veri miktarı)
        Xtype   % Karar değişkeni tipi % Unused
        LengthApp  % Her uygulamanın görev sayısını tutan vektör
    end
    
    % METHODS
    methods(Access = protected)
        %% ExtractParameters: Extract parameters from MegaGraph cell array
        function [nTask, Line, Column, S, T, W, Type, LengthApp] = ExtractParameters(obj)
            % Extract from obj.parameter (MegaGraph cell array)
            % MegaGraph format: {Task, Line, Column, S, T, W, Type, GraphLengths}
            
            if isempty(obj.parameter)
                error('MCMAOneTask: parameter (MegaGraph) is not provided');
            end
            
            param = obj.parameter;
            
            % Check if parameter is a cell array
            if ~iscell(param)
                error('MCMAOneTask: parameter must be a cell array (MegaGraph)');
            end
            
            % Extract parameters
            if length(param) < 7
                error('MCMAOneTask: MegaGraph must have at least 7 elements');
            end
            
            nTask = param{1};           % Task count (73)
            Line = param{2};            % NoC rows (9)
            Column = param{3};          % NoC columns (9)
            S = param{4};               % Source task indices
            T = param{5};               % Target task indices
            W = param{6};               % Edge weights
            Type = param{7};            % Type (usually 9)
            
            % GraphLengths is optional (8th element if provided)
            if length(param) >= 8
                LengthApp = param{8};
            else
                LengthApp = [];
            end
        end
    end
    
    methods
        %% Setting: Problem boyutunu tanımla
        function Setting(obj)
            [obj.nTask, obj.Line, obj.Column, obj.S, obj.T, obj.W, obj.Xtype, obj.LengthApp] = obj.ExtractParameters();
            
            if isempty(obj.M) 
                obj.M = 2;  % Amaç sayısı = 2 (Latency, Fault Tolerance)
            end
            
            if isempty(obj.D)
                obj.D = obj.nTask; % Karar değişkeni sayısı = görev sayısı
            end
            
            % Çözüm uzayı : Her görev [1, 81] arasında bir çekirdeğe atanmalı
            obj.lower = ones(1, obj.D);                             % Minimum: çekirdek 1
            obj.upper = (obj.Line * obj.Column) * ones(1, obj.D);   % Maksimum: çekirdek 81
            obj.encoding = obj.Xtype * ones(1, obj.D); 
        end
        
        %% Initialization: İlk popülasyonu oluştur
        function Population = Initialization(obj, N)
            if nargin < 2   % fonksiyon eksik parametre almış ise, gereksiz aslında
                N = obj.N;  % Popülasyon büyüklüğü = 100
            end
            
            Task = obj.nTask;
            CoresRow = obj.Line;
            CoresColumn = obj.Column;
            AppTaskLengths = obj.LengthApp;
            
            % İlk popülasyonu oluştur
            PopDec = InitializePopulation(N, Task, CoresRow, CoresColumn, AppTaskLengths);
            
            % Değerlendir
            Population = obj.Evaluation(PopDec);
        end
        
        %% CalObj: Her çözümü değerlendir
        function PopObj = CalObj(obj, PopDec)
            % Giriş: PopDec = N x 73 matrix (100 çözüm)
            % Çıkış: PopObj = N x 2 matrix (100 çözümün amaçları)

            % Parametreler
            Source = obj.S;  % Source görevler
            Target = obj.T;  % Target görevler
            Weight = obj.W;  % Edge ağırlıkları
            
            % NoC parametreleri
            nR = obj.Line;      % 9
            nC = obj.Column;    % 9
            
            % Çekirdek konumlarını hesapla
            [LN, CL] = ind2sub([nR nC], 1:nR*nC);
            Pos_Tab = [LN' CL']; % 81 x 2 (her çekirdeğin koordinatları)
            
            % Mesafe matrisini hesapla (Manhattan)
            Dist_Tab = pdist2(Pos_Tab, Pos_Tab, 'cityblock');
            
            % Cost hesaplamaları
            latency = CostCalculator.LatencyCost(PopDec, Dist_Tab, Source, Target, Weight);
            faultTolerance = CostCalculator.FaultToleranceCost(PopDec, nR, nC, Dist_Tab);
            
            PopObj = [latency, faultTolerance];
        end 
    end
end