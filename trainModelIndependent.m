% Trains a independent maximum-entropy model on a set of samples. This is a sub-function 
% invoked by trainModel() and typically should not be invoked directly.
%
% Usage:
%   model_out = trainModelExhaustive(model,samples)
%
% Input:
%   refer to trainModel.m
%
% Output:
%   model_out - trained ME model (normalized and with the entropy stored in model_out.H)
%
% Last update: Ori Maoz 30/06/2016
function model_out = trainModelIndependent(model,raster)


if (~strcmp(model.type,'indep'))
   error('This function can only train independent models'); 
end


% check if our input was a raster (of samples from a distribution) or an actual Boltzmann distribution
if (isstruct(raster))
    % we got a model (that we can compute the marginals directly on)
    use_exact_marginals = true;
    model_base = raster;
    ncells = model_base.ncells;
else
    use_exact_marginals = false;
    [ncells,nsamples] = size(raster);    
    if (nsamples < ncells)    
        warning('Input raster must be of the form (nsamples x ncells), are you sure you did not transpose it?');
    end
end



if (ncells ~= model.ncells)
    error('Number of cells in model does not match number of cells in data');
    
end


% compute the observables which are the mean firing rates. If we are preseted with experimental data,
% just take the empirical count. Otherwise we need to compute it explicitly on the base model
if (use_exact_marginals)
    npatterns = 2^ncells;
    unique_words = logical(de2bi(0:(npatterns-1)))';
    
    % get probabilities of the input patterns
    lp = mexLogProbability(unique_words,model_base);
    lp = lp - logsumexp(lp(:));
    empirical_probs = exp(lp);
    
    % get the observables - our goal is to fit the model marginals to these.
    firing_rates = mexGetMarginals(unique_words,model,empirical_probs)';

else    
    % mean firing rates
    firing_rates = mean(raster,2);
end


model_out = model;

% compute the partition function: z = (1-p1)(1-p2).....(1-pn)
logz = sum(log(1-firing_rates));
model_out.z = logz;


% now compute the entropy
% we need special handling for cells that never firing
nonzero_indices = firing_rates > 0;
model_out.H = sum(-firing_rates(nonzero_indices) .* log2(firing_rates(nonzero_indices)) - (1-firing_rates(nonzero_indices)) .* log2(1-firing_rates(nonzero_indices)));


% compute each of the factors
for i = 1:ncells
    
    if (firing_rates(i) > 0)    
        % first compute (1-p1)(1-p2)....(pi)(1-p_1+1)...(1-pn)
        vec = ones(ncells,1);
        vec(i) = 2*firing_rates(i);
        lambda = -(sum(log(vec - firing_rates)) - logz);
    else
       % this cell never fired, its Lagrange factor will be 0 
       % (the normal computation would result in infinity)
       lambda = 0;
        
    end
    model_out.factors(i) = lambda;
end

end