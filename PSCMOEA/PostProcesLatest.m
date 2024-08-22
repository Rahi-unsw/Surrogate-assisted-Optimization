%% General post process script for multiobjective problems
function PostProcesLatest
clear all;close all;clc;
filename = mfilename;
motherfolderpath = which(filename);
motherfolder = fileparts(motherfolderpath);
cd(motherfolder);
addpath(genpath(motherfolder));
FolderName = 'Data';
ProblemSet = 'MW'; % 'MW','LIRCMOP','DASCMOP','CF','FCP','Test', 'Combined'
% ProblemSet1 = 'Combined2';
DataType = 'Fig'; % 'Mat1', 'Mat2', 'Fig'
CompType = 'StateOfTheArt'; % 'StateOfTheArt', 'Variants', 'Illustration'

if strcmp(DataType,'Mat1')
    MoveData = strcat('ProcessedData',filesep,CompType,filesep,ProblemSet);
elseif strcmp(DataType,'Mat2')
    MoveData = strcat('ProcessedData',filesep,CompType,filesep,ProblemSet,filesep,'Data');
elseif strcmp(DataType,'Fig')
    MoveData = strcat('ProcessedData',filesep,CompType,filesep,ProblemSet,filesep,'Data');
    MoveFigure = strcat('ProcessedData',filesep,CompType,filesep,ProblemSet,filesep,'Figures');
end
DataFolder = strcat('ProcessedData',filesep,CompType,filesep,ProblemSet);

% 'weldedbeam', 'twobartruss', 'speedreducer', 'geartrain', 'discbrake', 'conceptualmarine', 'carsideimpact'
% Problems: 'MW1','MW2','MW3','MW4','MW5','MW6','MW7','MW8','MW9','MW10','MW11','MW12','MW13','MW14' - D = 10; M = 3 (MW4, MW8, MW14); M = 2 (MW1-MW3, MW5-MW7, MW9-MW13)
% 'CF1','CF2','CF3','CF4','CF5','CF6','CF7','CF8','CF9','CF10' - D = 10; M = 3 ('CF8','CF9','CF10'); M = 2 ('CF1','CF2','CF3','CF4','CF5','CF6','CF7')
% 'DASCMOP1','DASCMOP2','DASCMOP3','DASCMOP4','DASCMOP5','DASCMOP6','DASCMOP7','DASCMOP8','DASCMOP9' - D = 10; M = 3 ('DASCMOP7','DASCMOP8','DASCMOP9'); M = 2 ('DASCMOP1','DASCMOP2','DASCMOP3','DASCMOP4','DASCMOP5','DASCMOP6')
% 'LIRCMOP1','LIRCMOP2','LIRCMOP3','LIRCMOP4','LIRCMOP5','LIRCMOP6','LIRCMOP7','LIRCMOP8','LIRCMOP9','LIRCMOP10','LIRCMOP11','LIRCMOP12','LIRCMOP13','LIRCMOP14' - D = 10; M = 3 ('LIRCMOP13','LIRCMOP14'); M = 2 ('LIRCMOP1','LIRCMOP2','LIRCMOP3','LIRCMOP4','LIRCMOP5','LIRCMOP6','LIRCMOP7','LIRCMOP8','LIRCMOP9','LIRCMOP10','LIRCMOP11','LIRCMOP12')
% 'FCP1','FCP3','Test1','Test2'
% 'Test1','Test2','MW1','MW2','MW3','MW4','MW5','MW6','MW7','MW8','MW9','MW10','MW11','MW12','MW13','MW14','LIRCMOP1','LIRCMOP2','LIRCMOP3','LIRCMOP4','DASCMOP4','DASCMOP5','DASCMOP6','DASCMOP7','DASCMOP8','DASCMOP9'
% 'LIRCMOP5','LIRCMOP6','LIRCMOP7','LIRCMOP8','LIRCMOP9','LIRCMOP10','LIRCMOP11','LIRCMOP12','LIRCMOP13','LIRCMOP14','DASCMOP1','DASCMOP2','DASCMOP3'
param.allprobs = {'Test1','Test2','MW1','MW2','MW3','MW4','MW5','MW6','MW7','MW8','MW9','MW10','MW11','MW12','MW13','MW14','LIRCMOP1','LIRCMOP2','LIRCMOP3','LIRCMOP4','DASCMOP4','DASCMOP5','DASCMOP6','DASCMOP7','DASCMOP8','DASCMOP9'};
param.allruns = [31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31];% 31,31,31,31,31,31,31,31,31,31,31,31,31,31
param.algorithms = {'PSCMOEA'};% 'PSCMOEA','KMGSAEA','ASAMOEAD','KTS','MultiObjectiveEGO','SILE'
param.MFE_all = [500];
param.metric_all = {'IGD'}; % 'IGD', 'HV', 'IGDp', 'GD'

if strcmp(DataType,'Fig')
    MarkerEdgeColors = colormap('lines');
    linestyles = cellstr(char('-',':','-.','--','-',':','-.','--','-',':','-',':',...
        '-.','--','-',':','-.','--','-',':','-.','-',':','-.','--','-',':','-.','--'));
    Markers = {'o','+','*','x','s','d','^','v','>','<','p'};
end

%% Postprocessed data storing
if strcmp(DataType,'Mat1')
%     for mfe = 1:length(param.MFE_all)
%         param.MFE = param.MFE_all(mfe);
%         for i = 1:numel(param.allprobs)
%             param.prob_name = param.allprobs{i};
%             SplitStr = split(param.prob_name,"_");
%             prob = feval(param.prob_name);
%             [P] = PF(SplitStr{1},1000,prob.nf);
%             for j = 1:length(param.algorithms)
%                 for metric = 1:length(param.metric_all)
%                     param.metric = param.metric_all{metric};
%                     score_run = [];FF_FE = [];count = 0;
%                     for k = 1:param.allruns(i)
%                         if strcmp(FolderName,'Data')
%                             pth = strcat(cd,filesep,FolderName,filesep,param.algorithms{j},filesep,param.prob_name,filesep,'Trial-',num2str(k));
%                         end
%                         load(strcat(pth,filesep,param.prob_name,'_Archive.mat'));
%                         objective = Archive.muf(1:min(param.MFE,size(Archive.muf,1)),:); % param.MFE
%                         %                         if ~strcmp(ProblemSet,'LIRCMOP')
%                         [~,Archive.mug] = feval(param.prob_name,Archive.x);
%                         %                         end
%                         constraint = Archive.mug(1:min(param.MFE,size(Archive.mug,1)),:);
%                         violation = sum(max(constraint,0),2);
%                         NFE = 1:min(param.MFE,size(Archive.muf,1));
%                         idfeas = find(violation==0);
%                         if ~isempty(idfeas)
%                             f_feas = objective(idfeas,:);
%                             % Calculate performance metrics
%                             score = PerformMetric(param.metric,f_feas,P);
%                             % Score of each run
%                             score_run(k) = score;
%                             if ~exist(strcat(pwd,filesep,param.prob_name,'_',param.algorithms{j},'_NFE_',num2str(param.MFE),'_NFE_Statistics.mat'),'file')
%                                 % NFE to identify first feasible solution
%                                 firstfeasid = min(idfeas);
%                                 FF_FE(k) = NFE(firstfeasid);
%                             end
%                             count = count + 1;
%                         else
%                             score_run(k) = NaN;
%                             if ~exist(strcat(pwd,filesep,param.prob_name,'_',param.algorithms{j},'_NFE_',num2str(param.MFE),'_NFE_Statistics.mat'),'file')
%                                 FF_FE(k) = param.MFE; % max(NFE);
%                             end
%                         end
%                     end
%                     SC = count; % SR = 100*(count/param.allruns(i));
%                     save(strcat(param.prob_name,'_',param.algorithms{j},'_NFE_',num2str(param.MFE),'_',param.metric,'_RawScore_Statistics.mat'),'score_run');
%                     if ~exist(strcat(pwd,filesep,param.prob_name,'_',param.algorithms{j},'_NFE_',num2str(param.MFE),'_NFE_Statistics.mat'),'file')
%                         save(strcat(param.prob_name,'_',param.algorithms{j},'_NFE_',num2str(param.MFE),'_NFE_Statistics.mat'),'FF_FE');
%                         save(strcat(param.prob_name,'_',param.algorithms{j},'_NFE_',num2str(param.MFE),'_SuccessCount.mat'),'SC');
%                     end
%                 end
%             end
%         end
%     end


    % This part, I have to fire after storing the processed data and commenting that
    %% Finding the worst particular metric value of all runs of all algorithms compared for a particular problem
    for mfe = 1:length(param.MFE_all)
        param.MFE = param.MFE_all(mfe);
        for i = 1:numel(param.allprobs)
            param.prob_name = param.allprobs{i};
            for metric = 1:length(param.metric_all)
                param.metric = param.metric_all{metric};
                AllAlgoMetricStat = [];
                for j = 1:length(param.algorithms)
                    load(strcat(DataFolder,filesep,param.prob_name,'_',param.algorithms{j},'_NFE_',num2str(param.MFE),'_',param.metric,'_RawScore_Statistics.mat'),'score_run');
                    AllAlgoMetricStat = [AllAlgoMetricStat,score_run];
                end
                if strcmp(param.metric,'HV')
                    WorstMetricValue = max(min(AllAlgoMetricStat) - 1e-1,0); % eps = 1e-1
                else
                    WorstMetricValue = max(AllAlgoMetricStat) + 1e-1;
                end
                save(strcat(param.prob_name,'_NFE_',num2str(param.MFE),'_',param.metric,'_WorstValue.mat'),'WorstMetricValue');
            end
        end
    end




    %% Replace all the NAN metric value with the worst particular metric value
    for mfe = 1:length(param.MFE_all)
        param.MFE = param.MFE_all(mfe);
        for i = 1:numel(param.allprobs)
            param.prob_name = param.allprobs{i};
            for metric = 1:length(param.metric_all)
                param.metric = param.metric_all{metric};
                load(strcat(param.prob_name,'_NFE_',num2str(param.MFE),'_',param.metric,'_WorstValue.mat'),'WorstMetricValue');
                for j = 1:length(param.algorithms)
                    load(strcat(DataFolder,filesep,param.prob_name,'_',param.algorithms{j},'_NFE_',num2str(param.MFE),'_',param.metric,'_RawScore_Statistics.mat'),'score_run');
                    score_run(isnan(score_run)) = WorstMetricValue;
                    save(strcat(param.prob_name,'_',param.algorithms{j},'_NFE_',num2str(param.MFE),'_',param.metric,'_Score_Statistics.mat'),'score_run');
                end
            end
        end
    end
end



%% AllAlgo statistics for comparison
if strcmp(DataType,'Mat2')
    for mfe = 1:length(param.MFE_all)
        param.MFE = param.MFE_all(mfe);
        for metric = 1:length(param.metric_all)
            param.metric = param.metric_all{metric};
            for i = 1:numel(param.allprobs)
                param.prob_name = param.allprobs{i};
                for j = 1:length(param.algorithms)
                    load(strcat(DataFolder,filesep,param.prob_name,'_',param.algorithms{j},'_NFE_',num2str(param.MFE),'_',param.metric,'_Score_Statistics.mat'),'score_run');
                    [~,sid] = sort(score_run);
                    medianrun = sid((param.allruns(i)+1)/2);
                    each_stat(j,:) = [nanmin(score_run),nanmean(score_run),nanmedian(score_run),nanmax(score_run),nanstd(score_run),medianrun];
                    if ~exist(strcat(pwd,filesep,param.prob_name,'_NFE_',num2str(param.MFE),'_AllAlgo_NFEStat.mat'),'file')
                        load(strcat(DataFolder,filesep,param.prob_name,'_',param.algorithms{j},'_NFE_',num2str(param.MFE),'_NFE_Statistics.mat'),'FF_FE');
                        load(strcat(DataFolder,filesep,param.prob_name,'_',param.algorithms{j},'_NFE_',num2str(param.MFE),'_SuccessCount.mat'),'SC');
                        NFE_stat(j,:) = [nanmin(FF_FE),round(nanmean(FF_FE)),nanmedian(FF_FE),nanmax(FF_FE),round(nanstd(FF_FE)),SC];
                    end
                end
                Prob_stat = each_stat;
                save(strcat(param.prob_name,'_NFE_',num2str(param.MFE),'_',param.metric,'_AllAlgo_Statistics.mat'),'Prob_stat');
                if ~exist(strcat(pwd,filesep,param.prob_name,'_NFE_',num2str(param.MFE),'_AllAlgo_NFEStat.mat'),'file')
                    Prob_NFE = NFE_stat;
                    save(strcat(param.prob_name,'_NFE_',num2str(param.MFE),'_AllAlgo_NFEStat.mat'),'Prob_NFE');
                end
            end
        end
    end




    %% Statistical significance test
    for mfe = 1:length(param.MFE_all)
        param.MFE = param.MFE_all(mfe);
        for metric = 1:length(param.metric_all)
            param.metric = param.metric_all{metric};
            for i = 1:numel(param.allprobs)
                param.prob_name = param.allprobs{i};
                for j = 1:length(param.algorithms)
                    load(strcat(DataFolder,filesep,param.prob_name,'_',param.algorithms{j},'_NFE_',num2str(param.MFE),'_',param.metric,'_RawScore_Statistics.mat'),'score_run');
                    sig_data_each(:,j) = score_run';
                    if ~exist(strcat(pwd,filesep,'Overall_Significance_Stat_NFE_','NFE_',num2str(param.MFE),'.mat'),'file')
                        load(strcat(DataFolder,filesep,param.prob_name,'_',param.algorithms{j},'_NFE_',num2str(param.MFE),'_NFE_Statistics.mat'),'FF_FE');
                        sig_data_each1(:,j) = FF_FE';
                    end
                end
                if strcmp(param.metric,'HV')
                    sig_data{i} = -sig_data_each;
                else
                    sig_data{i} = sig_data_each;
                end
                if ~exist(strcat(pwd,filesep,'Overall_Significance_Stat_NFE_','NFE_',num2str(param.MFE),'.mat'),'file')
                    if strcmp(ProblemSet,'DASCMOP')
                        sig_data_each1(sig_data_each1<50) = 1;
                    end
                    sig_data1{i} = sig_data_each1;
                end
            end
            [OverallSigStat,PairwiseEachAlgo] = Wilcoxon(sig_data);
            save(strcat('Pairwise_Significance_Each_Algo_','NFE_',num2str(param.MFE),'_',param.metric,'.mat'),'PairwiseEachAlgo');
            save(strcat('Overall_Significance_Stat_','NFE_',num2str(param.MFE),'_',param.metric,'.mat'),'OverallSigStat');
            if ~exist(strcat(pwd,filesep,'Overall_Significance_Stat_NFE_','NFE_',num2str(param.MFE),'.mat'),'file')
                [OverallSigStatNFE,PairwiseEachAlgoNFE] = Wilcoxon(sig_data1);
                save(strcat('Pairwise_Significance_Each_Algo_NFE_','NFE_',num2str(param.MFE),'.mat'),'PairwiseEachAlgoNFE');
                save(strcat('Overall_Significance_Stat_NFE_','NFE_',num2str(param.MFE),'.mat'),'OverallSigStatNFE');
            end
        end
    end
end



if strcmp(DataType,'Fig')
    %     %% Performance profile
    %     for mfe = 1:length(param.MFE_all)
    %         param.MFE = param.MFE_all(mfe);
    %         for metric = 1:length(param.metric_all)
    %             param.metric = param.metric_all{metric};
    %             for i = 1:numel(param.allprobs)
    %                 param.prob_name = param.allprobs{i};
    %                 load(strcat(MoveData,filesep,param.prob_name,'_NFE_',num2str(param.MFE),'_',param.metric,'_AllAlgo_Statistics.mat'),'Prob_stat');
    %                 if ~exist(strcat(pwd,filesep,'Performance_Profile_Meanrun_',CompType,'_',ProblemSet,'_NFE_',num2str(param.MFE),'.jpg'),'file')
    %                     load(strcat(MoveData,filesep,param.prob_name,'_NFE_',num2str(param.MFE),'_AllAlgo_NFEStat.mat'),'Prob_NFE');
    %                     MeanAllNFE(:,i) = Prob_NFE(:,2);
    %                 end
    %                 if strcmp(param.metric,'HV')
    %                     Mean_all(:,i) = -Prob_stat(:,2);
    %                 else
    %                     Mean_all(:,i) = Prob_stat(:,2);
    %                 end
    %             end
    %             PerformanceProfile(Mean_all,param.algorithms);
    %             filename = strcat('Performance_Profile_Meanrun_',CompType,'_',ProblemSet,'_NFE_',num2str(param.MFE),'_',param.metric);
    %             saveas(gcf,filename);
    %             print(filename,'-depsc2','-r300');
    %             saveas(gcf,strcat(filename,'.jpg'));
    %             if ~exist(strcat(pwd,filesep,'Performance_Profile_Meanrun_',CompType,'_',ProblemSet,'_NFE_',num2str(param.MFE),'.jpg'),'file')
    %                 MeanAllNFE(MeanAllNFE<10) = 1;
    %                 PerformanceProfile(MeanAllNFE,param.algorithms);
    %                 filename1 = strcat('Performance_Profile_Meanrun_',CompType,'_',ProblemSet,'_NFE_',num2str(param.MFE));
    %                 saveas(gcf,filename1);
    %                 print(filename1,'-depsc2','-r300');
    %                 saveas(gcf,strcat(filename1,'.jpg'));
    %             end
    %         end
    %     end





    % %% Convergence plots
    % for mfe = 1:length(param.MFE_all)
    %     param.MFE = param.MFE_all(mfe);
    %     for metric = 1:length(param.metric_all)
    %         param.metric = param.metric_all{metric};
    %         for i = 1:numel(param.allprobs)
    %             param.prob_name = param.allprobs{i};
    %             SplitStr = split(param.prob_name,"_");
    %             prob = feval(param.prob_name);
    %             [P] = PF(SplitStr{1},1000,prob.nf);
    %             for j = 1:length(param.algorithms)
    %                 load(strcat(DataFolder,filesep,param.prob_name,'_',param.algorithms{j},'_NFE_',num2str(param.MFE),'_',param.metric,'_Score_Statistics.mat'),'score_run');
    %                 [~,sid] = sort(score_run);
    %                 medianrun = sid((param.allruns(i)+1)/2);
    %                 if strcmp(FolderName,'Data')
    %                     pth = strcat(cd,filesep,'Data',filesep,param.algorithms{j},filesep,param.prob_name,filesep,'Trial-',num2str(medianrun));
    %                 end
    %                 load(strcat(pth,filesep,param.prob_name,'_Archive.mat'));
    %                 NFE = 1:min(param.MFE,size(Archive.muf,1));
    %                 objective = Archive.muf(NFE,:);
    %                 [~,Archive.mug] = feval(param.prob_name,Archive.x);
    %                 constraint = Archive.mug(NFE,:);
    %                 violation = sum(max(constraint,0),2);
    %                 idfeas = find(violation==0);
    %                 idinfeas = setdiff(NFE,idfeas,'stable');
    %                 % Calculating convergence data
    %                 if ~isempty(idfeas)
    %                     ConData = [];
    %                     if isempty(idinfeas)
    %                         init = 1; % min(size(f,1),11*prob.nf - 1);
    %                         count = init;
    %                         while count <= length(NFE)
    %                             ff = objective(1:count,:);
    %                             SC = PerformMetric(param.metric,ff,P);
    %                             if count > init
    %                                 SC = ConValid(param.metric,ConData(:,2),SC);
    %                             end
    %                             ConData = [ConData;[count SC]];
    %                             count = count + 5;
    %                         end
    %                         SC = PerformMetric(param.metric,objective,P);
    %                         if ConData(end,1) < param.MFE
    %                             ConData = [ConData;[param.MFE SC]];
    %                         end
    %                         FConData = interp1(ConData(:,1),ConData(:,2),1:param.MFE);
    %                     else
    %                         FFeasID = idfeas(1);
    %                         ff = objective(FFeasID,:);
    %                         SC = PerformMetric(param.metric,ff,P);
    %                         ConData = [ConData;[FFeasID SC]];
    %                         if length(idfeas) > 1
    %                             for count = 2:length(idfeas)
    %                                 ff = objective(idfeas(1:count),:);
    %                                 SC = PerformMetric(param.metric,ff,P);
    %                                 SC = ConValid(param.metric,ConData(:,2),SC);
    %                                 ConData = [ConData;[idfeas(count) SC]];
    %                             end
    %                         end
    %                         if ConData(end,1) < param.MFE
    %                             ConData = [ConData;[param.MFE SC]];
    %                         end
    %                         FConData = interp1(ConData(:,1),ConData(:,2),1:param.MFE);
    %                     end
    %                     plot(1:param.MFE,FConData,'DisplayName',param.algorithms{j},'color',MarkerEdgeColors(j,:),'Marker',Markers{j},'MarkerSize',7,'MarkerIndices',1:5:param.MFE,'linestyle',linestyles{j},'LineWidth',2);hold on;
    %                     %                 xlim([0 NFEvec(end)]);
    %                 else
    %                     MetricVal = NaN(1,param.MFE);
    %                     plot(1:param.MFE,MetricVal,'DisplayName',param.algorithms{j},'color',MarkerEdgeColors(j,:),'Marker',Markers{j},'MarkerSize',7,'MarkerIndices',1:5:param.MFE,'linestyle',linestyles{j},'LineWidth',2);hold on;
    %                 end
    %             end
    %             xlabel('Evaluation count');ylabel(strcat('Median-run',{' '},param.metric,{' '},'value'));
    %             set(gca,'YScale', 'log' ,'FontSize',14,'FontName','Arial');
    %             legend('Location','northeast','FontSize',12);legend('show');grid on;
    %             hold off;
    %             filename = strcat(param.prob_name,'_NFE_',num2str(param.MFE),'_',param.metric,'_AllAlgo_ConvergencePlot');
    %             saveas(gcf,filename);
    %             print(filename,'-depsc2','-r300');
    %             saveas(gcf,strcat(filename,'.jpg'));
    %         end
    %     end
    % end




            %% SearchSwitching plots
        %     'MW','MW','MW','MW','MW','MW','MW','MW','MW','MW','MW','MW','MW','MW'
        %     ,'Test','Test','LIRCMOP','LIRCMOP','LIRCMOP','LIRCMOP','DASCMOP','DASCMOP','DASCMOP','DASCMOP','DASCMOP','DASCMOP'
            ProblemSet = {'Test','Test','MW','MW','MW','MW','MW','MW','MW','MW','MW','MW','MW','MW','MW','MW',...
                'LIRCMOP','LIRCMOP','LIRCMOP','LIRCMOP','DASCMOP','DASCMOP','DASCMOP','DASCMOP','DASCMOP','DASCMOP'};
            for mfe = 1:length(param.MFE_all)
                param.MFE = param.MFE_all(mfe);
                for metric = 1:length(param.metric_all)
                    param.metric = param.metric_all{metric};
                    for i = 1:numel(param.allprobs)
                        param.prob_name = param.allprobs{i};
                        DataFolder = strcat('ProcessedData',filesep,CompType,filesep,ProblemSet{i});
                        for j = 1:length(param.algorithms)
                            load(strcat(DataFolder,filesep,param.prob_name,'_',param.algorithms{j},'_NFE_',num2str(param.MFE),'_',param.metric,'_Score_Statistics.mat'),'score_run');
                            [~,sid] = sort(score_run);
                            medianrun = sid((param.allruns(i)+1)/2);
                            if strcmp(FolderName,'Data')
                                pth = strcat(cd,filesep,'Data',filesep,param.algorithms{j},filesep,param.prob_name,filesep,'Trial-',num2str(medianrun));
                            end
                            load(strcat(pth,filesep,param.prob_name,'_Archive.mat'));
                            UnconsSearch = sum(Archive.SearchFlag == 1);
                        end
                        AllUnconsSearch(i) = UnconsSearch;
                    end
                    % Plot
                    bar(AllUnconsSearch);hold on;
                    xtickangle(45);
                    hAx = gca;
                    hAx.XTickLabelRotation = 45;
                    hAx.XTick = 1:length(param.allprobs);
                    hAx.XTickLabel = param.allprobs;
                    ylim([0 390]);
                    ylabel('Iterations after initialization');
                    set(gca,'FontSize',18,'FontName','Arial');
%                     f = gcf;f.InnerPosition = [300 300 1500 700]
                    hold off;
    %                 filename = strcat('TestLIRDASSeries_SearchSwitching');
                    filename = strcat('TestMWLIRDASSeries_SearchSwitching');
                    saveas(gcf,filename);
                    print(filename,'-depsc2','-r300');
                    saveas(gcf,strcat(filename,'.jpg'));
                end
            end


%     %% Figures to show the converged solutions in objective space
%     for mfe = 1:length(param.MFE_all)
%         param.MFE = param.MFE_all(mfe);
%         for metric = 1:length(param.metric_all)
%             param.metric = param.metric_all{metric};
%             for i = 1:numel(param.allprobs)
%                 param.prob_name = param.allprobs{i};
%                 prob = feval(param.prob_name);
% 
%                 % Visualization data generation
%                 SplitStr = split(param.prob_name,"_");
%                 [P] = PF(SplitStr{1},1000,prob.nf);
%                 [RR] = PFUnCon(SplitStr{1},1000,prob.nf);
%                 R = FR(param.prob_name,prob.nf);
%                 if ~isempty(R)
%                     xx = R{1};yy = R{2};zz = R{3};
%                 else
%                     xx = [];yy = []; zz = [];
%                 end
% 
%                 for j = 1:length(param.algorithms)
%                     load(strcat(DataFolder,filesep,param.prob_name,'_',param.algorithms{j},'_NFE_',num2str(param.MFE),'_',param.metric,'_Score_Statistics.mat'),'score_run');
%                     [~,sid] = sort(score_run);
%                     medianrun = sid((param.allruns(i)+1)/2);
%                     for k = 1:param.allruns(i)
%                         medianrun = 26;
%                         %                     medianrun = 9; % 9 for MW11
%                         pth = strcat(cd,filesep,'Data',filesep,param.algorithms{j},filesep,param.prob_name,filesep,'Trial-',num2str(medianrun));
%                         load(strcat(pth,filesep,param.prob_name,'_Archive.mat'));
%                         NFE = 1:min(param.MFE,size(Archive.muf,1));
%                         objective = Archive.muf(NFE,:);
%                         %                     if ~strcmp(ProblemSet,'LIRCMOP')
%                         [~,Archive.mug] = feval(param.prob_name,Archive.x);
%                         %                     end
%                         constraint = Archive.mug(NFE,:);
%                         violation = sum(max(constraint,0),2);
%                         idfeas = find(violation==0);
%                         FeasF = objective(idfeas,:);
%                         idinfeas = setdiff(NFE,idfeas,'stable');
%                         InFeasF = objective(idinfeas,:);
% 
%                         if prob.nf < 3
%                             %                         % Visualization
%                             % %                         rng(101,'twister');
%                             % %                         pop = Initialize_pop(prob,100);
%                             % %                         [muf,mug] = feval(param.prob_name,pop.x);
%                             %                         figure;
%                             %                         s1 = surf(xx,yy,zz);hold on;
%                             %                         view(0,90);
%                             %                         s1.EdgeColor = 'none';s1.FaceAlpha = 0.7;
%                             %                         p1 = plot(RR(:,1),RR(:,2),'k.','LineWidth',2);
%                             %                         p2 = plot(P(:,1),P(:,2),'r.','LineWidth',2);
%                             %                         p3 = plot(muf(:,1),muf(:,2),'k.');
%                             %                         legend([s1,p1,p2,p3],{'FR','UPF','CPF','Initial sol.'},'Orientation','horizontal');
%                             %                         box on;ax = gca;ax.BoxStyle = 'full';
%                             %                         xlabel('f_1');ylabel('f_2');
%                             %                         %                 title(strcat(param.prob_name,'\_',param.algorithms{j}));
%                             %                         set(gca,'fontsize',18);
%                             %                         hold off;
%                             %                         filename = strcat(param.prob_name,'_case3');
%                             %                         saveas(gcf,filename);
%                             %                         print(filename,'-depsc2','-r300');
%                             %                         saveas(gcf,strcat(filename,'.jpg'));
% 
% %                             figure;
% %                             s1 = surf(xx,yy,zz);hold on;
% %                             view(0,90);
% %                             s1.EdgeColor = 'none';s1.FaceAlpha = 0.7;s1.DisplayName = 'FR';
% %                             plot(RR(:,1),RR(:,2),'k.','DisplayName','UPF');
% %                             plot(P(:,1),P(:,2),'r.','DisplayName','CPF');
% %                             plot(InFeasF(:,1),InFeasF(:,2),'k.','MarkerSize',10,'DisplayName','Infeas.');
% %                             plot(FeasF(:,1),FeasF(:,2),'b.','MarkerSize',10,'DisplayName','Feas.');
% %                             box on;ax = gca;ax.BoxStyle = 'full';legend('show');
% %                             xlabel('f_1');ylabel('f_2');
% %                             %                 title(strcat(param.prob_name,'\_',param.algorithms{j}));
% %                             set(gca,'fontsize',16);
% %                             hold off;
% % %                             filename = strcat(param.prob_name,'_',param.algorithms{j},'_AllSamples');
% % %                             saveas(gcf,filename);
% % %                             print(filename,'-depsc2','-r300');
% % %                             saveas(gcf,strcat(filename,'.jpg'));
% 
%                             figure;
%                             s1 = surf(xx,yy,zz);hold on;
%                             view(0,90);
%                             s1.EdgeColor = 'none';s1.FaceAlpha = 0.7;s1.DisplayName = 'FR';
%                             plot(RR(:,1),RR(:,2),'k.','DisplayName','UPF');
%                             plot(P(:,1),P(:,2),'r.','DisplayName','CPF');
%                             plot(InFeasF(:,1),InFeasF(:,2),'k.','MarkerSize',10,'DisplayName','Infeas.');
%                             plot(FeasF(:,1),FeasF(:,2),'b.','MarkerSize',10,'DisplayName','Feas.');
%                             box on;ax = gca;ax.BoxStyle = 'full';legend('show');
%                             xlabel('f_1');ylabel('f_2');
%                             %                 title(strcat(param.prob_name,'\_',param.algorithms{j}));
%                             set(gca,'fontsize',16);
%                             PlotLimits(param.prob_name,P);
%                             hold off;
%                             filename = strcat(param.prob_name,'_',param.algorithms{j},'_AllSamples_Zoomed');
%                             saveas(gcf,filename);
%                             print(filename,'-depsc2','-r300');
%                             saveas(gcf,strcat(filename,'.jpg'));
% 
%                         else
%                             %                     figure;
%                             %                     p1 = plot3(RR(:,1),RR(:,2),RR(:,3),'k.');hold on;
%                             %                     p2 = plot3(P(:,1),P(:,2),P(:,3),'r.');
%                             %                     legend([p1,p2],{'Unconstrained PF','Constrained PF'});
%                             %                     box on;ax = gca;ax.BoxStyle = 'full';
%                             %                     xlabel('f_1');ylabel('f_2');
%                             %                     %                 title(strcat(param.prob_name,'\_',param.algorithms{j}));
%                             %                     set(gca,'fontsize',16);
%                             %                     hold off;
%                             %                     filename = strcat(param.prob_name,'_properties');
%                             %                     saveas(gcf,filename);
%                             %                     print(filename,'-depsc2','-r300');
% 
%                             %                         figure;
%                             %                         %                     s1 = surf(xx,yy,zz);hold on;
%                             %                         %                     view(0,90);
%                             %                         %                     s1.EdgeColor = 'none';s1.FaceAlpha = 0.7;
%                             %                         plot3(RR(:,1),RR(:,2),RR(:,3),'k.');hold on;
%                             %                         plot3(P(:,1),P(:,2),P(:,3),'r.');view(47,19);
%                             %                         plot3(InFeasF(:,1),InFeasF(:,2),InFeasF(:,3),'k.','MarkerSize',10);
%                             %                         plot3(FeasF(:,1),FeasF(:,2),FeasF(:,3),'b.','MarkerSize',10);
%                             %                         box on;ax = gca;ax.BoxStyle = 'full';
%                             %                         xlabel('f_1');ylabel('f_2');zlabel('f_3');
%                             %                         %                 title(strcat(param.prob_name,'\_',param.algorithms{j}));
%                             %                         set(gca,'fontsize',16);
%                             %                         hold off;
%                             %                         filename = strcat(param.prob_name,'_',param.algorithms{j},'_AllSamples');
%                             %                         saveas(gcf,filename);
%                             %                         print(filename,'-depsc2','-r300');
%                             %                         saveas(gcf,strcat(filename,'.jpg'));
%                             %                         %
%                             figure;
%                             %                     s1 = surf(xx,yy,zz);hold on;
%                             %                     view(0,90);
%                             %                     s1.EdgeColor = 'none';s1.FaceAlpha = 0.7;
%                             plot3(RR(:,1),RR(:,2),RR(:,3),'k.');hold on;
%                             plot3(P(:,1),P(:,2),P(:,3),'r.');view(47,19);
%                             plot3(InFeasF(:,1),InFeasF(:,2),InFeasF(:,3),'k.','MarkerSize',10);
%                             plot3(FeasF(:,1),FeasF(:,2),FeasF(:,3),'b.','MarkerSize',10);
%                             box on;ax = gca;ax.BoxStyle = 'full';
%                             xlabel('f_1');ylabel('f_2');zlabel('f_3');
%                             %                 title(strcat(param.prob_name,'\_',param.algorithms{j}));
%                             set(gca,'fontsize',16);
%                             PlotLimits(param.prob_name,P);
%                             hold off;
%                             filename = strcat(param.prob_name,'_',param.algorithms{j},'_AllSamples_Zoomed');
%                             saveas(gcf,filename);
%                             print(filename,'-depsc2','-r300');
%                             saveas(gcf,strcat(filename,'.jpg'));
%                         end
%                     end
%                 end
%             end
%         end
%     end
end

%% Move files
if strcmp(DataType,'Mat1') || strcmp(DataType,'Mat2')
    movefile('*.mat',MoveData);
else
    movefile('*.fig',MoveFigure);
    movefile('*.eps',MoveFigure);
    movefile('*.jpg',MoveFigure);
end
return


%% Initializing a population
function pop = Initialize_pop(prob,no_solutions)
pop.x = repmat(prob.bounds(:,1)',no_solutions,1)+repmat((prob.bounds(:,2)-prob.bounds(:,1))',no_solutions,1).*lhsdesign(no_solutions,prob.nx);
return



%% Plot limits
function PlotLimits(Prob,P)
switch Prob
    case 'MW1'
        xlim([0 1]);ylim([0 2]);
    case 'MW2'
        xlim([0 1]);ylim([0 2]);
    case 'MW3'
        xlim([0 1]);ylim([0 2]);
    case 'MW4'
        xlim([0 1.5]);ylim([0 1.5]);zlim([0 1.5]);
    case 'MW5'
        xlim([0 2]);ylim([0 2]);
    case 'MW6'
        xlim([0 1.5]);ylim([0 2]);
    case 'MW7'
        xlim([0 2]);ylim([0 2]);
    case 'MW8'
        xlim([0 1.5]);ylim([0 1.5]);zlim([0 1.5]);
    case 'MW9'
        xlim([0 2]);ylim([0 2]);
    case 'MW10'
        xlim([0 1]);ylim([0 2]);
    case 'MW11'
        xlim([0 3]);ylim([0 3]);
    case 'MW12'
        xlim([0 2]);ylim([0 2]);
    case 'MW13'
        xlim([0 2]);ylim([0 5]);
    case 'MW14'
        xlim([0 1.5]);ylim([0 1.5]);zlim([0 6]);
    case 'Test1'
        xlim([0 10]);ylim([0 10]);
    case 'Test2'
        xlim([0 10]);ylim([0 10]);
    case 'DASCMOP1'
        xlim([0 2.5]);ylim([0 2.5]);
    case 'DASCMOP2'
        xlim([0 2.5]);ylim([0 2.5]);
    case 'DASCMOP3'
        xlim([0 2.5]);ylim([0 2.5]);
    case 'LIRCMOP5'
        xlim([0 5]);ylim([0 5]);
    case 'LIRCMOP6'
        xlim([0 5]);ylim([0 5]);
    case 'LIRCMOP9'
        xlim([0 5]);ylim([0 5]);
    case 'LIRCMOP10'
        xlim([0 5]);ylim([0 5]);
    case 'LIRCMOP11'
        xlim([0 5]);ylim([0 5]);
    otherwise
        if size(P,2) < 3
            xlim([min(P(:,1)) max(P(:,1))]);ylim([min(P(:,2)) max(P(:,2))]);
        else
            xlim([min(P(:,1)) max(P(:,1))]);ylim([min(P(:,2)) max(P(:,2))]);zlim([min(P(:,3)) max(P(:,3))]);
        end
end
return


%% Select performance metric
function score = PerformMetric(metric,f,P)
switch metric
    case 'IGD'
        score = IGD(f,P);
    case 'GD'
        score = GD(f,P);
    case 'HV'
        %         score = HV(f,P);
        score = HV_Revised(f,P);
    case 'IGDp'
        %         score = IGD_Plus(f,P);
        score = IGDp(f,P);
end
return


%% Select Convergence validation
function score = ConValid(metric,data,score)
if strcmp(metric,'IGD') || strcmp(metric,'IGD_plus')
    if score <= min(data)
        score = score;
    else
        score = min(data);
    end
elseif strcmp(metric,'HV')
    if score >= max(data)
        score = score;
    else
        score = max(data);
    end
end
return