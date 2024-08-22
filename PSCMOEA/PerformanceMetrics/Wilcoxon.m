%% Significance test
function [OverallSigStat,PairwiseEachAlgo] = Wilcoxon(Data)
% --------------------Input structure----------------------
%   Data: a cell matrix, cell size = no of problems
%   Each cell containt M*N matrix where,
%   M = no of run; N = no of algorithms
% --------------------Output structure----------------------
%   Final: a M*N matrix where M = no of algorithms and
%   N = [Win Tie Loose]; It contains the overall pairwise 
%   significance comparison of any algorithm with all others
%   across all problems.
%   PairwiseEachAlgo:  = Significance stat of each algorithm
%   (cell size = no of algorithms) with others (rowwise) across
%   each problem (columnwise)
OverallSigStat = [];
PairwiseEachAlgo = cell(1,size(Data{1},2));
for m = 1:size(Data{1},2)
    a = [];
    for n = 1:numel(Data)
        ranksum_test_input = Data{n};
        base_algo = ranksum_test_input(:,m);
        for o = 1:size(ranksum_test_input,2)
            input = [base_algo,ranksum_test_input(:,o)];
            [~,sigResult] = wilcoxenranksum_significance_test(input);
            SigTestResult{n,1}{o,1} = o;
            SigTestResult{n,1}{o,2} = sigResult;
        end
        tmp = cell2mat(SigTestResult{n}) ;
        a = [a tmp(:,2)];
        win(n) = length(find(tmp(:,2)==1));
        loose(n) = length(find(tmp(:,2)==-1));
        tie(n) = length(find(tmp(:,2)==0));
    end
    PairwiseEachAlgo{m} = a;
    tmp = [win;tie;loose];
    final = [sum(win) sum(tie)-numel(Data) sum(loose)];
    OverallSigStat = [OverallSigStat;final];
end
end


%% Wilcoxen Ranksum Significance Test
function [hall,sigResult] = wilcoxenranksum_significance_test(ranksum_test_input)
% Please cite, 'Wilcoxon, Frank, S. K. Katti, and Roberta A. Wilcox. 
% Critical values and probability levels for the Wilcoxon rank sum test 
% and the Wilcoxon signed rank test. Pearl River, NY: American Cyanamid Company, 1963'.
hall = [];
if ~isempty(ranksum_test_input)
    igdalgo(:,1) = ranksum_test_input(:,1);
    igdalgo(:,2) = ranksum_test_input(:,2);
    worstvalue = max([igdalgo(:,1)',igdalgo(:,2)']) + 1e-1; % eps = 1e-1
    igdalgo(isnan(igdalgo)) = worstvalue;
    try
        [~,h1] = ranksum(igdalgo(:,1),igdalgo(:,2),'alpha',0.05,'tail','left');
        [~,h2] = ranksum(igdalgo(:,1),igdalgo(:,2),'alpha',0.05,'tail','right');
    catch
        if sum(isnan(igdalgo(:,1)))==size(ranksum_test_input,1) && sum(isnan(igdalgo(:,2)))~=size(ranksum_test_input,1)
            h1 = 0;
            h2 = 1;
        elseif sum(isnan(igdalgo(:,1)))~=size(ranksum_test_input,1) && sum(isnan(igdalgo(:,2)))==size(ranksum_test_input,1)
            h1 = 1;
            h2 = 0;
        elseif sum(isnan(igdalgo(:,1)))==size(ranksum_test_input,1) && sum(isnan(igdalgo(:,2)))==size(ranksum_test_input,1)
            h1 = 0;
            h2 = 0;
        end
    end
    hall = [hall,h1 h2];
    if h1 == 0 && h2 == 0
        sigResult = 0;
    elseif h1 == 0 && h2 == 1
        sigResult = -1;
    elseif h1 == 1 && h2 == 0
        sigResult = 1;
    else
        sigResult = ' ';
    end
else
    sigResult = ' ';
end
end