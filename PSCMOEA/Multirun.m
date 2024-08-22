function Multirun
%% Parallel Computing
clc;clear all;close all;warning off;
filename = mfilename;
motherfolderpath = which(filename);
motherfolder = fileparts(motherfolderpath);
cd(motherfolder);
addpath(genpath(motherfolder));
Params;path = cd;count = 1;
isparallel = 0; % 1 means parfor; 0 means no parallel

for i = 1:numel(param.allprobs)
    param.prob_name = param.allprobs{i}; % param.prob = param.probs{i};
    prob = feval(param.prob_name);
    if prob.nf < 3
        param.spacing = 99;
    elseif prob.nf == 3
        param.spacing = 12;
    end
    for j = 1:length(param.algorithms)
        param.algo = param.algorithms{j};
        for k = 32:param.allruns(i)
            disp(strcat(param.algo,', Problem Name:- ',param.prob_name,', Trial:- ',num2str(k)));
            pth = strcat(cd,filesep,'ExpData',filesep,param.algo,filesep,param.prob_name,filesep,'Trial-',num2str(k));
            mkdir(pth);
            datapath = [pth,filesep,param.prob_name,'_Archive.mat'];
            if ~exist(datapath,'file')
                cd(pth);
                param.seed = 100+k;
                param.MFE = param.MFE_all(i);
                save('Parameters.mat', 'param');
                path1{count} = pth;
                count = count+1;
                cd(path);
            end
        end
    end
end
cd(path);


clear param;
%% parfor: utilize multiple core of PC.
if ~isempty(path1)
    if isparallel == 1
        numcores = feature('numcores');
        parpool(numcores);
        parfor i = 1:length(path1)
            cd(path1{i});
            disp(strcat('Running -> ',path1{i}));
            tic;
            Driver(path1{i});
            toc;
        end
        delete(gcp);
    else
        for i = 1:length(path1)
            cd(path1{i});
            disp(strcat('Running -> ',path1{i}));
            tic;
            Driver(path1{i});
            toc;
        end
    end
end
cd(path);
return