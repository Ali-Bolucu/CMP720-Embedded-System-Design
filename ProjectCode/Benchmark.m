%% Benchmark.m


%% Configurations Setup
% ======================================================
addpath(fileparts(mfilename('fullpath')));

scriptDir = fileparts(mfilename('fullpath'));
if ~isempty(scriptDir)
    addpath(scriptDir);
    cd(scriptDir);
end

Npop   = 100;
OutDir = fullfile(scriptDir, 'results', 'figures');
if ~exist(OutDir, 'dir'); mkdir(OutDir); end
%%


%% Create MegaGraph
% ======================================================
MegaGraph = ApplicationToMegaGraph();
%%
    
%% Generate Base Random Pool
global BASE_RANDOM_POOL;
NumTasks = MegaGraph{1};
NumRouters = MegaGraph{2} * MegaGraph{3};
BASE_RANDOM_POOL = zeros(Npop, NumTasks);
for i = 1:Npop
    BASE_RANDOM_POOL(i, :) = randperm(NumRouters, NumTasks);
end

USE_SPECTRAL_FLAG = true;
CONVERGENCE_TRACE = struct('latency', [], 'thermal', []);
%% Run PlatEMO
% ======================================================
[Dec, Obj, Con] = platemo( ...
    'N',         Npop,                  ...  % population size
    'maxFE',     100000,                ...  % maximum function evaluations
    'problem',   @NocMappingProblem,    ...  % Problem function
    'algorithm', @NSGAII,               ...  % 
    'parameter', MegaGraph,             ...  % 
    'save',      1                      ...  % 
);
%% 


%% Visualize Results
% ======================================================
% Extract best solution from Pareto front
BestDec = Dec(1, :);      % First solution (typically best)
BestObj = Obj(1, :);      % Corresponding objectives

% Generate and save figures
RunLabel = datestr(now, 'yyyymmdd_HHMMSS');
VisualizeResults(Obj, BestDec, BestObj, MegaGraph, OutDir, RunLabel);

fprintf('[Benchmark] Results saved to: %s\n', OutDir);
%%

cd(scriptDir);