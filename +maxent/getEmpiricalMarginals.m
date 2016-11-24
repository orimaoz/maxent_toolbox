% Returns the empirical marginals of a ME model for a set of samples
%
% Usage:
%   marginals = getEmpiricalMarginals(samples,model)
%
%
% The function accepts an array of samples and a model and returns the empirical average of the model's moment
% functions over the set of samples.
%
% Arguments (mandatory):
%   samples - Set of samples to train the model on, in the format (ncells x nsamples). 
%   model   - Maximum entropy model as returned by trainModel(). This model is only used to choose the appropriate set
%             of marginals, the actual model parameters (lagrange multipliers in the model.factors) are irrelevant.
%
% Output:
%   marginals - an (1xnfactors) vector of marginals.
%
%
% Last update: Ori Maoz 10/07/2016
function marginals = getEmpiricalMarginals(samples, model)

if (nargin<2)
    error('Usage: marginals = getEmpiricalMarginals(samples,model)');
end


% delegate to the MEX implementation
marginals = mexEmpiricalMarginals(samples,model);

end