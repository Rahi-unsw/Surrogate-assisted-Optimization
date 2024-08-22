function SigTableGenerator
clc;clear all;close all;warning off;
filename = mfilename;
motherfolderpath = which(filename);
motherfolder = fileparts(motherfolderpath);
cd(motherfolder);
addpath(genpath(motherfolder));
ProblemSet = 'DASCMOP'; % 'MW','LIRCMOP','DASCMOP','CF'
CompType = 'StateOfTheArt'; % 'StateOfTheArt', 'Variants'
pathmat = strcat(pwd,filesep,'ProcessedData',filesep,CompType,filesep,ProblemSet);

Metric = {'IGD','IGDp','HV','NFE'};FEMax = 500;% FolderName = 'MOBO';,'KMGSAEA','ASAMOEAD','KTS','MultiObjectiveEGO','SILE'
ColLabels = {'SASSCMOEA','KMGSAEA','ASAMOEAD','KTS','MultiObjectiveEGO','SILE'};
% 'weldedbeam', 'twobartruss', 'speedreducer', 'geartrain', 'discbrake', 'conceptualmarine', 'carsideimpact'
% 'MW1','MW2','MW3','MW4','MW5','MW6','MW7','MW8','MW9','MW10','MW11','MW12','MW13','MW14'
Problemnames = {'DASCMOP'};
SetRows = 1; % Calculate
Alignment = 'c'; % 'c', 'l'
FileName = strcat('Sig_stat_AllMetric','.tex');
TableLabel = 'tab:Performance comparison';

cd(pathmat)
if ~exist('Tables','dir')
    mkdir Tables;
end
cd Tables;

%% Generate the table
fid = fopen(FileName, 'w');
fprintf(fid,'\\begin{table*}[!ht]\\scriptsize\n');%\\scriptsize
fprintf(fid,'\\centering\n');
fprintf(fid,'\\caption{Statistical significance test comparison on %s series.}\n',ProblemSet);
fprintf(fid,'\\label{%s}\n',TableLabel);
fprintf(fid,'\\tabcolsep 1.3mm\n');
fprintf(fid,'\\renewcommand{\\arraystretch}{1.7}\n');
fprintf(fid,'\\begin{center}\n\\begin{tabular}{');

fprintf(fid, '%c%c%c', Alignment, Alignment, Alignment);
for i = 1:numel(ColLabels) %width
    fprintf(fid, '%c', Alignment);
end
fprintf(fid, '}\n');
fprintf(fid, '\\hline\n');

if(~isempty(ColLabels))
    fprintf(fid, '\\multirow{2}{*}{\\textbf{Metric}} & \\multirow{2}{*}{\\textbf{Problems}} & \\multirow{2}{*}{\\textbf{M}} & \\multicolumn{%d}{c}{\\textbf{Algorithms}}\\\\\\cline{4-%d} & & \n',numel(ColLabels),(numel(ColLabels)+3));
    for w = 1:numel((ColLabels))-1
        fprintf(fid, ' & \\multicolumn{1}{c}{\\textbf{%s}}', ColLabels{w});
    end
    fprintf(fid, ' & \\multicolumn{1}{c}{\\textbf{%s}}\\\\\\hline\n', ColLabels{end});
end
fprintf(fid, '\r\n');

% Precess data
for mm = 1:numel(Metric)
    fprintf(fid, '\\multirow{%d}{*}{\\textbf{%s}} ', SetRows, Metric{mm});
    for p = 1:numel(Problemnames)
        P = Problemnames{p};
        if ~isempty(regexp(P,'MW', 'once')) || ~isempty(regexp(P,'CF', 'once')) || ~isempty(regexp(P,'LIRCMOP', 'once')) || ~isempty(regexp(P,'DASCMOP', 'once'))
            RowLabels = [2]; % [2 3 5 7]
        end
        fprintf(fid, '& \\multirow{%d}{*}{\\textbf{%s}} ', length(RowLabels), P);
        count = 1;
        for m = RowLabels
            if ~isempty(regexp(P,'MW', 'once')) || ~isempty(regexp(P,'CF', 'once')) || ~isempty(regexp(P,'LIRCMOP', 'once')) || ~isempty(regexp(P,'DASCMOP', 'once'))
                if strcmpi(Metric{mm},'NFE')
                    load(strcat(pathmat,filesep,'Data',filesep,'Overall_Significance_Stat_NFE_NFE_',num2str(FEMax),'.mat'));
                    win1 = OverallSigStatNFE(:,1)'; Tie1 = OverallSigStatNFE(:,2)'; Loss1 = OverallSigStatNFE(:,3)';
                else
                    load(strcat(pathmat,filesep,'Data',filesep,'Overall_Significance_Stat_NFE_',num2str(FEMax),'_',Metric{mm},'.mat'));
                end
            end
            win = OverallSigStat(:,1)'; Tie = OverallSigStat(:,2)'; Loss = OverallSigStat(:,3)';
            
            % Create the problem cells
            if strcmpi(Metric{mm},'NFE')
                if count == 1
                    fprintf(fid, '& %d ', m);
                else
                    fprintf(fid, '& & %d ', m);
                end
                for w = 1:numel(ColLabels)
                    fprintf(fid, '& %d/%d/%d ', win1(w), Tie1(w), Loss1(w));
                end
            else
                if count == 1
                    fprintf(fid, '& %d ', m);
                else
                    fprintf(fid, '& & %d ', m);
                end
                for w = 1:numel(ColLabels)
                    fprintf(fid, '& %d/%d/%d ', win(w), Tie(w), Loss(w));
                end
            end
            fprintf(fid, '\\\\\n');
            count = count + 1;
        end
    end
    fprintf(fid,'\\hline\n\n');
end
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{center}\n');
fprintf(fid, '\\end{table*}\r\n');
fclose(fid);
cd(motherfolder);
end