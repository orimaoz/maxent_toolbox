% Exhaustively calculate the marginals of a maximum entropy model.
% This function goes over all the 2^ncells states, and so it is not recommended to run it on very large models.
%
% Usage:
%   marginals = getMarginals(model)
%
% The function accepts a model and an array of samples and returns the log probability for each sample 
% according to the model. If the model has not been normalize, the results will also be un-normalized.
%
% Arguments (mandatory):
%   model   - Maximum entropy model as returned by trainModel()
%
% Output:
%   marginals - an (1xnfactors) vector of marginals.
%
%
% Last update: Ori Maoz 30/06/2016
function marginals = getMarginals(model)

if (model.ncells > 30)
    error('Your model is big. This is a heavy computation. Are you sure you want to do this?');
end

% Now delegate to the MEX implementation
marginals = mexExhaustiveMarginals(model);

end