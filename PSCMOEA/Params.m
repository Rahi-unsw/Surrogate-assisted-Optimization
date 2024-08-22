%% Parameters to run the algorithm
% Offspring generation scheme and parameters
param.prob_crossover = 0.9;
param.prob_mutation = 0.1;
param.distribution_crossover = 10;
param.distribution_mutation = 20;

% Running parameters
% 'weldedbeam', 'twobartruss', 'speedreducer', 'geartrain', 'discbrake', 'conceptualmarine', 'carsideimpact'
% Problems: 'MW1','MW2','MW3','MW4','MW5','MW6','MW7','MW8','MW9','MW10','MW11','MW12','MW13','MW14' - D = 10; M = 3 (MW4, MW8, MW14); M = 2 (MW1-MW3, MW5-MW7, MW9-MW13)
% 'CF1','CF2','CF3','CF4','CF5','CF6','CF7','CF8','CF9','CF10' - D = 10; M = 3 ('CF8','CF9','CF10'); M = 2 ('CF1','CF2','CF3','CF4','CF5','CF6','CF7')
% 'DASCMOP1','DASCMOP2','DASCMOP3','DASCMOP4','DASCMOP5','DASCMOP6','DASCMOP7','DASCMOP8','DASCMOP9' - D = 10; M = 3 ('DASCMOP7','DASCMOP8','DASCMOP9'); M = 2 ('DASCMOP1','DASCMOP2','DASCMOP3','DASCMOP4','DASCMOP5','DASCMOP6')
% 'LIRCMOP1','LIRCMOP2','LIRCMOP3','LIRCMOP4','LIRCMOP5','LIRCMOP6','LIRCMOP7','LIRCMOP8','LIRCMOP9','LIRCMOP10','LIRCMOP11','LIRCMOP12','LIRCMOP13','LIRCMOP14' - D = 10; M = 3 ('LIRCMOP13','LIRCMOP14'); M = 2 ('LIRCMOP1','LIRCMOP2','LIRCMOP3','LIRCMOP4','LIRCMOP5','LIRCMOP6','LIRCMOP7','LIRCMOP8','LIRCMOP9','LIRCMOP10','LIRCMOP11','LIRCMOP12')
param.allprobs = {'Test1'}; % Specify the problems to fire
param.allruns = [32]; % Number of runs for each problem
param.MFE_all = [500]; % Specify NFE for each problem
% param.strategies = [1]; % Specify the variants you want to fire
param.algorithms = {'PSCMOEA'}; % Specify the variants you want to fire
param.visualization = 1; % 1 means ON and 0 means OFF

% Algorithmic parameters
param.EA_popsize = 100;
param.EA_generations = 100;
param.x_threshold = 1e-4;
param.x_neighbour = [];
param.infeas_ratio = 0.2;

