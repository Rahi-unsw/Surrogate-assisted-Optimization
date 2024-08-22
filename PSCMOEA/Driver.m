%% To fire all variants
function Driver(path)
% Load the parameters;
load('Parameters.mat');
if strcmpi(param.algo,'SASSCMOEA_test')
    SASSCMOEA_test(path);           % Variant 1;
elseif strcmpi(param.algo,'PSCMOEA')
    PSCMOEA(path);
end
return

