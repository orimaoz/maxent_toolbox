% Samples from a maximum entropy model
%
% Usage:
%   samples = generateSamples(model, nsamples)
%   samples = generateSamples(model, nsamples, Name,Value,...)
%
% The function accepts a model and the number of required samples.
% according to the model. If the model has not been normalize, the results will also be un-normalized.
%
% Arguments (mandatory):
%   model    - maximum entropy model as returned by trainModel()
%   nsamples - how many samples to generate
%   separation - how many samples to jump each time (default = take every sample, 2 = ignore every other sample,
%                3 = take one sample, drop the two next samples etc etc)
%
% Optional arguments (in the form of Name,Value pairs):
%   x0   - starting state for MCMC walk
%
% Output:
%   logprobs - an (1xnsamples) vector of log probabilities in natural logarithm.
%
%
% Last update: Ori Maoz 30/06/2016
function samples = generateSamples(model, nsamples,varargin)

p = inputParser;
addOptional(p,'x0',0);          % starting state
addOptional(p,'burnin',10000);  % number of burn-in samples
addOptional(p,'separation',1);  % number of burn-in samples
p.parse(varargin{:});
x0 = p.Results.x0;
burnin = p.Results.burnin;
separation = p.Results.separation;


if (numel(x0)>1)
    if (numel(x0) ~= model.ncells)
        error('starting state did not have the same number of dimensions as model');
    end
    
    % TODO make a copy of it in memory? Currently it might change in-place in some scenarios
end



% delegate to the MEX implementation
params.burnin = burnin;
params.separation = separation;
samples = gibbsSampler(x0,nsamples,model,params);


end