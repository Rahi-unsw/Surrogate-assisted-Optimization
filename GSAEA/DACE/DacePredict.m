function [y, sigma] = DacePredict(X, dmodel)
% if nargin ~= 2
% 	error('dace_predict requires 2 input arguments')
% end
samples = size(X,1);
if samples > 1
    [f, ssqr] = predictor(X, dmodel);
else
    [f, ~, ssqr] = predictor(X, dmodel);
end
if size(f,1) == size(X,1)
    y = f;
else
    y = f';
end
sigma = sqrt(abs(ssqr));
return
