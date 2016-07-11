% Samples from a maximum entropy model
%
% Usage:
%   samples = generateSamples(model, nsamples)
%   samples = generateSamples(model, nsamples, Name,Value,...)
%
% The function accepts a model and the number of samples to be generated from it.
% The samples are generated via a Metropolis-Hastings sampling scheme. In order to obtain samples which
% represent the target disribution, default settings include dropping the first 10000 generated samples ("burn-in") and skipping 
% most of the samples generated along the way ("separation"). These parameters can be controlled by the user.
%
% Arguments (mandatory):
%   model    - maximum entropy model as returned by trainModel()
%   nsamples - how many samples to generate
%
% Optional arguments (in the form of Name,Value pairs):
%   x0      - starting state for MCMC walk. This can also be "0" signifying the all-zero (0000...) pattern.
%   burnin  - how many samples to drop before returning results (default: 10000).
%   separation - how many bit flips to perform before taking each sample (one for every dimension of the input).
%
% Output:
%   samples - an (ncells x nsamples) matrix of samples from the model.
%
%
% Last update: Ori Maoz 30/06/2016
function samples = generateSamples(model, nsamples,varargin)

DEFAULT_BURNIN = 10000;

p = inputParser;
addOptional(p,'x0',0);          % starting state
addOptional(p,'burnin',DEFAULT_BURNIN);  % number of burn-in samples
addOptional(p,'separation',model.ncells);  % number of burn-in samples
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
samples = mexGibbsSampler(model,nsamples,x0,params);


end