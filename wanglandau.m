% Estimates the partition function and entropy of Boltzmann distributions on binary inputs using the Wang-Landau method.
% Parameters: 
%   model:      Boltzmann distribution to work on. At the time of writing this, the implementation supported only Ising
%               model and MERP model.
%
%   Optional arguments (in string, value pairs).
%   binsize:    Bin size for the energy histogram. Default if not inputed is 0.01. 
%   depth:      Accuracy parameter for the simulation (integer). Higher values mean a higher accuracy and longer runtime.
%               The final accuracy is in the order of exp(2^-(depth-1))
%   separation: Number of samples to skip for every sample obtained in the MCMC random walk. A larger value decorrelates
%               the samples and provides more accurate results, but incurs a longer run-time.
%   analytic:   If this value is set as 'true' then the contents of the energy bins are computed analytically instead of
%               using MCMC methods. This is mainly used to evaluate the accuracy of the MCMC process and will only work
%               on small (up to ~25 dimension) distributions.
%
% Output:
%   adds to the model the fields 'z' which is the log partition function and 'H' which is the entropy (in bits).
function model_out = wanglandau(model,varargin);

DEFAULT_BIN_SIZE = 0.01;
DEFAULT_DEPTH = 20;
DEFAULT_SEPARATION = 1;
DEFAULT_ANALYTIC = false;

if nargin==0
    error('usage: wanglandau(model,[optional arguments])')
    
end


% parse optional arguments
p = inputParser;
addOptional(p,'binsize',DEFAULT_BIN_SIZE,@isnumeric); % energy bin size
addOptional(p,'depth',DEFAULT_DEPTH,@isnumeric); % depth of MCMC interations
addOptional(p,'separation',DEFAULT_SEPARATION,@isnumeric); % separation between samples
addOptional(p,'analytic',DEFAULT_ANALYTIC,@islogical);    % exhaustive computation
p.parse(varargin{:});
binsize = p.Results.binsize;
depth = p.Results.depth;
separation = p.Results.separation;
analytic = p.Results.analytic;


% remove any current value of the partition function, we will estimate it by ourselves
model.z = 0;

if (analytic)
    % exhaustive computation of the energy histogram
    [e, g] = energyDensity(model,binsize);
else
    % Use MCMC on the energies
    t=tic;
    [e, g] = WLS_estimateDensity(model,binsize,depth,separation);
    toc(t)
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
model_out.H = H;


end