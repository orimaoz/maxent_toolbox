% Returns the empirical marginals of a ME model for a set of samples
%
% Usage:
%   marginals = getEmpiricalMarginals(model, samples)
%
%
% The function accepts a model and an array of samples and returns the empirical average of the model's moment
% functions over the set of samples.
%
% Arguments (mandatory):
%   model   - Maximum entropy model as returned by trainModel()
%   samples - Set of samples to train the model on, in the format (ncells x nsamples). 
%
% Output:
%   marginals - an (1xnfactors) vector of marginals.
%
%
% Last update: Ori Maoz 30/06/2016
function marginals = getEmpiricalMarginals(model, samples)

% delegate to the MEX implementation
marginals = mexGetMarginals(samples,model);

end