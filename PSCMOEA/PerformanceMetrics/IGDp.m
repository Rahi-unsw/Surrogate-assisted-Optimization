function Score = IGDp(PopObj,PF)
fmin   = min(PF,[],1);fmax   = max(PF,[],1);
PopObj = (PopObj-fmin)./(fmax-fmin);
PF = (PF-fmin)./(fmax-fmin);
% <metric> <min>
% Inverted generational distance Plus

if size(PopObj,2) ~= size(PF,2)
    Score = nan;
else
    [Nr,M] = size(PF);
    [N,~]  = size(PopObj);
    delta  = zeros(Nr,1);
    for i = 1 : Nr
        delta(i) = min(sqrt(sum(max(PopObj - repmat(PF(i,:),N,1),zeros(N,M)).^2,2)));
    end
    Score = mean(delta);
end
end