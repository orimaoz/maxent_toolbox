% Returns the marginals of a maximum entropy model, by either exhaustive computation (for small models) or MCMC 
% sampling (for big models).
%
% Usage:
%   marginals = getMarginals(model)
%   marginals = getMarginals(model,Name,Value,...)
%
% Arguments (mandatory):
%   model   - Maximum entropy model as returned by trainModel()
%
% Optional arguments (in the form of Name,Value pairs):
%   force_exhaustive   - force exhaustive computation. This could incur a very long runtime for big models.
%   force_mcmc         - force MCMC (sampling-based) computation, which returns only an approximate result.
%   nsamples           - number of empirical sampled used to compute the marginals for the MCMC-based method
%
% Output:
%   marginals - an (1xnfactors) vector of marginals.
%
% Last update: Ori Maoz 05/07/2016
function marginals = getMarginals(model,varargin)

DEFAULT_NSAMPLES = 100000; 

if (nargin<1)
    error('Usage: marginals = getMarginals(model,<Name,Value>,...)');
end


% check if the user wants to force any type of solver
p = inputParser;
addOptional(p,'force_exhaustive',false,@islogical);   % force exhaustive solver
addOptional(p,'force_mcmc',false,@islogical);         % force MCMC solver
addOptional(p,'nsamples',DEFAULT_NSAMPLES,@isnumeric);           % samples to use for MCMC marginal estimation
p.KeepUnmatched = true;
p.parse(varargin{:});


% should we get the marginals exhaustively or using MCMC?
solve_exhaustive = (model.ncells < 30);
if (p.Results.force_mcmc)
    solve_exhaustive = false;
end
if (p.Results.force_exhaustive)
    if (model.ncells > 30)
        warning('Your model is big. This is a heavy computation. Are you sure you want to do this?');        
    end
    solve_exhaustive = true;
end


if (solve_exhaustive)
    % Exhaustive computation - delegate to the MEX implementation
    marginals = mexExhaustiveMarginals(model);
else
    % sample from the model
    samples = maxent.generateSamples(model,p.Results.nsamples);
    
    % return empirical marginals on the sample set
    marginals = maxent.getEmpiricalMarginals(samples,model);
end



end