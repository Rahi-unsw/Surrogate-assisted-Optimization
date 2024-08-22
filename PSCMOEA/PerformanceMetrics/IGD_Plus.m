function Score = IGD_Plus(PopObj,PF)
fmin   = min(PF,[],1);fmax   = max(PF,[],1);
PopObj = (PopObj-fmin)./(fmax-fmin);
PF = (PF-fmin)./(fmax-fmin);
% <metric> <min>
% Inverted generational distance Plus
for i = 1:size(PF,1)
    ref = repmat(PF(i,:),size(PopObj,1),1);
    tmp = max((PopObj-ref),0);
    ff(i,:) = sqrt(sum(tmp.^2,2));
end
Distance_Plus = min(ff,[],2);
Score = mean(Distance_Plus);
end