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
%   x0      - starting state for MCMC walk. Default: the all-zero (0000...) pattern.
%   burnin  - how many samples to drop before returning results (default: 10000).
%   separation - how many bit flips to perform before taking each sample (default: one for every dimension of the input).
%   fix_indices - indices to elements that should be fixed during sampling (i.e. remain unchanged)
%
% Output:
%   samples - an (ncells x nsamples) matrix of samples from the model.
%
%
% Last update: Ori Maoz 09/11/2017
function samples = generateSamples(model, nsamples,varargin)

DEFAULT_BURNIN = 10000;

if nargin<2
    error('Usage: generateSamples(model, nsamples, ...)');
end

p = inputParser;
addOptional(p,'x0',0);          % starting state
addOptional(p,'burnin',DEFAULT_BURNIN);  % number of burn-in samples
addOptional(p,'separation',model.ncells);  % number of burn-in samples
addOptional(p,'randseed',nan);  % random seed
addOptional(p,'fix_indices',[]);  % indices to elements that should be fixed during sampling (i.e. remain unchanged)
p.parse(varargin{:});
x0 = p.Results.x0;
burnin = p.Results.burnin;
separation = p.Results.separation;
randseed = p.Results.randseed;
fix_indices = p.Results.fix_indices;


if (numel(x0)>1)
    if (numel(x0) ~= model.ncells)
        error('starting state did not have the same number of dimensions as model');
    end    
end

if (((max(fix_indices)) > model.ncells) | (min(fix_indices) < 1))
    error('fixed indices should be in the range of 1...#cells');
end


% change the "fixed" elements into a list of "to change" elements that the internal implementation expects.
% the "-1" is because the C code expects indices to be zero-based
params.variable_indices = uint32(setdiff(1:model.ncells,fix_indices) - 1);

if (numel(params.variable_indices) == 0)
    error('at least some of the indices should remain unfixed');
end


% delegate to the MEX implementation
if ~isnan(randseed)
    params.randseed = randseed;
end
params.burnin = burnin;
params.separation = separation;
samples = mexGibbsSampler(model,nsamples,x0,params);


end