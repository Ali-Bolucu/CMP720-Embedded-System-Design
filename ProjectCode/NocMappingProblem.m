%% NocMappingProblem.m
% Many-Core to NoC Task Mapping Problem
% 
% Optimizes task-to-router assignments on a 2D mesh NoC architecture.
% Minimizes: (1) Communication Cost, (2) Load Imbalance

classdef NocMappingProblem < PROBLEM
    
    properties(Access = private)
        nTasks              % Total tasks (e.g., 73 for 4 applications)
        Line                % NoC mesh rows (e.g., 9)
        Column              % NoC mesh columns (e.g., 9)
        S                   % Source task indices (edges)
        T                   % Target task indices (edges)
        W                   % Edge weights (communication volume)
        GraphLengths        % Cumulative task counts: [16 28 53 73]
        Xtype = 7           % Encoding type: 7 = permutation
        TamApp = []         % Application module sizes (from ParameterSet)
        useSpectral = true  % true: spectral seeding, false: random init
    end
    
    methods
        function Setting(obj)
            % Unpack parameters via PlatEMO's ParameterSet
            [obj.nTasks, obj.Line, obj.Column, obj.S, obj.T, obj.W, obj.Xtype, obj.TamApp] = ...
                obj.ParameterSet(73, 9, 9, [], [], [], 7, []);
            
            % Configure decision and objective space
            if isempty(obj.D); obj.D = obj.nTasks; end
            if isempty(obj.M); obj.M = 2; end
            
            % Variable bounds (baseline compatible: 0 to 1000)
            obj.lower = zeros(1, obj.D);
            obj.upper = 1000 * ones(1, obj.D);
            
            % Encoding type: 7 = permutation (baseline compatible)
            obj.encoding = obj.Xtype * ones(1, obj.D);
            
        fprintf('[NoCMappingProblem] nTasks=%d NoC=%dx%d Xtype=%d edges=%d TamApp=[%s]\n', ...
                obj.nTasks, obj.Line, obj.Column, obj.Xtype, length(obj.S), ...
                mat2str(obj.TamApp));
        end
        
        function Population = Initialization(obj, N)
            % Generate initial population using spatial locality-aware method
            if nargin < 2; N = obj.N; end
            
            numberofTasks = obj.nTasks;
            nocSize = obj.Line * obj.Column;
            dimensionApp = obj.TamApp;
            
            % Check for global override (used by CompareConvergence script)
            global USE_SPECTRAL_FLAG;
            if isempty(USE_SPECTRAL_FLAG)
                flag = obj.useSpectral;
            else
                flag = USE_SPECTRAL_FLAG;
            end
            
            % Use population creator (English version of baseline initialization)
            PopDec = PopulationCreator(N, numberofTasks, nocSize, dimensionApp, flag, obj.S, obj.T, obj.W);
            
            Population = obj.Evaluation(PopDec);
        end
        
        function PopObj = CalObj(obj, PopDec)
            % Calculate two objectives using baseline method (MCMACustoV4)
            
            % Generate NoC core coordinates
            nR = obj.Line;
            nC = obj.Column;
            [LN, CL] = ind2sub([nR nC], 1:nR*nC);
            Pos_Tab = [LN' CL'];
            
            % Calculate Manhattan distance matrix
            Dist_Tab = pdist2(Pos_Tab, Pos_Tab, 'cityblock');
            
            % Evaluate population using NoCObjectiveEvaluator (Latency and Thermal Hotspot)
            PopObj = NocObjectiveEvaluator(PopDec, Dist_Tab, nR*nC, ...
                                           obj.S, obj.T, obj.W, nR, nC);
        end
        
    end
    % ===== END OF CLASS =====
end
