% checks whether the model is within (threshold) standard deviations of the empirical marginals
function result = testMarginals(model, data, threshold)

nsamples = size(data,2);

% we assume that getEmpiricalMarginals and getMarginals have already been verified to work fine...
empirical_marginals = maxent.getEmpiricalMarginals(data,model);
model_marginals = maxent.getMarginals(model);

% find the expected standard error of each marginal
[phat, marginals_pci] = binofit(round(empirical_marginals * nsamples),nsamples,0.32);
empirical_marginals_std = max(abs(repmat(phat',[1 2]) - marginals_pci),[],2)';    

% normalize the errors to units of standard deviations (sample noise)
marginal_error = abs(model_marginals - empirical_marginals) ./ empirical_marginals_std;

result = max(marginal_error) < threshold;

end