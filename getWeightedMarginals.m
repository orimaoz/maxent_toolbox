% Returns the marginals of a ME model for a set of samples weighted by a set of probabilities
%
% Usage:
%   marginals = getWeightedMarginals(model, samples, probabilities)
%
%
% The function accepts a model and an array of samples and returns the empirical average of the model's moment
% functions over the set of samples.
%
% Arguments (mandatory):
%   model   - Maximum entropy model as returned by trainModel()
%   samples - Set of samples to train the model on, in the format (ncells x nsamples). 
%   probabilities - the probability used to reweight each sample
%
% Output:
%   marginals - an (1xnfactors) vector of marginals.
%
%
% Last update: Ori Maoz 30/06/2016
function marginals = getWeightedMarginals(model, samples,probabilities)

[ncells,nsamples] = size(samples);

if (numel(probabilities) ~= nsamples)
    error('number of samples must be equal to the number of probabilities');
end

% delegate to the MEX implementation
marginals = mexEmpiricalMarginals(samples,model,probabilities);

end