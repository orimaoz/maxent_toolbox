% Estimates the partition function and entropy of Boltzmann distributions on binary inputs using the Wang-Landau method.
%
% Usage:
%   model_out = wangLandau(model)
%   model_out = wangLandau(model, Name,Value,...)
%
% Arguments (mandatory):
%   model:      Maximum entropy model as returned by trainModel();
%
% Optional arguments (in the form of Name,Value pairs):
%   binsize:    Bin size for the energy histogram (Default: 0.01).
%   depth:      Accuracy parameter for the simulation (integer). Higher values mean a higher accuracy and longer runtime.
%               The final accuracy is in the order of exp(2^-(depth-1))
%   separation: Number of samples to skip for every sample obtained in the MCMC random walk. A larger value decorrelates
%               the samples and provides more accurate results, but incurs a longer run-time.
%   savefile         - will constantly save the state in this file, and try to resume from it if it already exists.
%   save_delay       - delay between saves (in seconds).
%
% Output:
% model_out - Model structure with two additional fields appended to it:
% model_out.z - Log partition function of the model.
% model_out.entropy -Entropy of the model (in bits).
%
% Last update: Ori Maoz, July 2017
%
function model_out = wangLandau(model,varargin)

DEFAULT_BIN_SIZE = 0.01;
DEFAULT_DEPTH = 20;
DEFAULT_SEPARATION = 1;
DEFAULT_ANALYTIC = false;

% number of seconds after we save our current status to a file
DEFAULT_TIME_BETWEEN_SAVES = 60;


if nargin==0
    error('usage: wanglandau(model,[optional arguments])')
    
end


% parse optional arguments
p = inputParser;
addOptional(p,'binsize',DEFAULT_BIN_SIZE,@isnumeric); % energy bin size
addOptional(p,'depth',DEFAULT_DEPTH,@isnumeric); % depth of MCMC interations
addOptional(p,'separation',DEFAULT_SEPARATION,@isnumeric); % separation between samples
addOptional(p,'analytic',DEFAULT_ANALYTIC,@islogical);    % exhaustive computation
addOptional(p,'savefile','');                 % save file for re-entrant code
addOptional(p,'save_delay',DEFAULT_TIME_BETWEEN_SAVES);    % time between re-entry file saves

p.parse(varargin{:});
binsize = p.Results.binsize;
depth = p.Results.depth;
separation = p.Results.separation;
analytic = p.Results.analytic;
savefile = p.Results.savefile;
save_delay = p.Results.save_delay;

% remove any current value of the partition function, we will estimate it by ourselves
model.z = 0;

if (analytic)
    % exhaustive computation of the energy histogram
    [e, g] = maxent.energyDensity(model,binsize);
else
    % Use MCMC on the energies
    [e, g] = maxent.WLS_estimateDensity(model,binsize,depth,separation,savefile,save_delay);
end

% since g is not normalized we can scale it down to save numerical problems.
% g it in log-scale so scaling it down is equivalent to subtraction.
% we have to take care to do it only in the non-zero entries of the histogram.
minval = min(g);
g = g - minval;

% g is an unnormalized log of energy state densities. Un-log it and normalize to 2^N:
g = exp(g);
g = g / sum(g);
g = g * 2^model.ncells(1);

% now we can estimate the partition function
Z = sum(g .* exp(-e));
logZ = log(Z);

% use the partition function to estimate the entropy
H = sum(e/log(2) .* g .* exp(-e))/sum(g.*exp(-e)) + log2(sum(g.*exp(-e)));

% return the partition function and the entropy in the model
model_out = model;
model_out.z = -logZ;
model_out.entropy = H;


end