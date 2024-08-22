function [y, sigma] = DacePredict(X, dmodel)
% if nargin ~= 2
% 	error('dace_predict requires 2 input arguments')
% end
samples = size(X,1);
[f, ssqr] = predictor(X, dmodel);
if size(f,1) == size(X,1)
    y = f;
else
    y = f';
end
sigma = sqrt(max(0,ssqr));
return
