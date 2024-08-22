%% Parameters to run the algorithm
% Offspring generation scheme and parameters
param.prob_crossover = 0.9;
param.prob_mutation = 0.1;
param.distribution_crossover = 10;
param.distribution_mutation = 20;

% Problem parameters
% 'TR1','TR2','TR3' :- TR 2 objectives problems
% 'MaF1_M2_D6','MaF2_M2_D6','MaF3_M2_D6','MaF4_M2_D6','MaF5_M2_D6','MaF7_M2_D6','MaF10_M2_D6','MaF11_M2_D7','MaF12_M2_D6' :- MaF 2 Objective problems
% 'MaF1_M3_D9','MaF2_M3_D9','MaF3_M3_D9','MaF4_M3_D9','MaF5_M3_D9','MaF6_M3_D9','MaF7_M3_D9','MaF8_M3_D2','MaF9_M3_D2','MaF10_M3_D9','MaF11_M3_D10','MaF12_M3_D9','MaF13_M3_D5' :- MaF 3 Objective problems
% 'MaF1_M5_D15','MaF2_M5_D15','MaF3_M5_D15','MaF4_M5_D15','MaF5_M5_D15','MaF6_M5_D15','MaF7_M5_D15','MaF8_M5_D2','MaF9_M5_D2','MaF10_M5_D15','MaF11_M5_D16','MaF12_M5_D15','MaF13_M5_D5' :- MaF 5 Objective problems
% 'MaF1_M7_D21','MaF2_M7_D21','MaF3_M7_D21','MaF4_M7_D21','MaF5_M7_D21','MaF6_M7_D21','MaF7_M7_D21','MaF8_M7_D2','MaF9_M7_D2','MaF10_M7_D21','MaF11_M7_D22','MaF12_M7_D21','MaF13_M7_D5' :- MaF 7 Objective problems
% 'DTLZ1_M2_D6','DTLZ2_M2_D6','DTLZ3_M2_D6','DTLZ4_M2_D6','DTLZ5_M2_D6','DTLZ6_M2_D6' :- DTLZ 2 Objective problems
% 'DTLZ1_M3_D9','DTLZ2_M3_D9','DTLZ3_M3_D9','DTLZ4_M3_D9','DTLZ5_M3_D9','DTLZ6_M3_D9' :- DTLZ 3 Objective problems
% 'DTLZ1_M5_D15','DTLZ2_M5_D15','DTLZ3_M5_D15','DTLZ4_M5_D15','DTLZ5_M5_D15','DTLZ6_M5_D15' :- DTLZ 5 Objective problems
% 'DTLZ1_M7_D21','DTLZ2_M7_D21','DTLZ3_M7_D21','DTLZ4_M7_D21','DTLZ5_M7_D21','DTLZ6_M7_D21' :- DTLZ 7 Objective problems
% 'ZDT1_M2_D6','ZDT2_M2_D6','ZDT3_M2_D6','ZDT4_M2_D6','ZDT6_M2_D6' :- ZDT 2 Objective problems

% Running parameters
param.allprobs = {'MaF1_M2_D6','MaF2_M2_D6','MaF3_M2_D6'}; % Specify the problems to fire
param.allruns = [2,2,2]; % Number of runs for each problem
param.MFE_all = [300,300,300]; % Specify NFE for each problem

% Algorithmic parameters
param.EA_popsize = 100;
param.EA_generations = 100;
param.iterate = 5;
param.x_threshold = 1e-4;
param.x_neighbour = [];