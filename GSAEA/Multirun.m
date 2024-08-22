function Multirun
%% Parallel Computing
clc;clear all;close all;warning off;
filename = mfilename;
motherfolderpath = which(filename);
motherfolder = fileparts(motherfolderpath);
cd(motherfolder);
addpath(genpath(motherfolder));
Params;path = cd;count = 1;

for i = 1:numel(param.allprobs)
    param.prob_name = param.allprobs{i}; % param.prob = param.probs{i};
    prob = feval(param.prob_name);
    if prob.nf < 3
        param.spacing = 99;
    elseif prob.nf == 3
        param.spacing = 12;
    elseif prob.nf == 5
        param.spacing = 4;
    elseif prob.nf == 7
        param.spacing = 3;
    end
    for k = 1:param.allruns(i)
        disp(strcat('GSAEA: Problem Name:- ',param.prob_name,', Trial:- ',num2str(k)));
        pth = strcat(cd,filesep,'Data',filesep,'GSAEA',filesep,param.prob_name,filesep,'Trial-',num2str(k));
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
cd(path);

clear param;
% parfor: utilize multiple core of PC.
numcores = feature('numcores');
parpool(numcores);
parfor i = 1:length(path1)
    cd(path1{i});
    disp(strcat('Running -> ',path1{i}));
    tic;
    GSAEA(path1{i});
    toc;
end
delete(gcp);
cd(path);
return