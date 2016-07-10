% Example code for the maximum entropy toolkit
% Ori Maoz, July 2016:  orimaoz@gmail.com,


%% part 1: working with small distributions of neurons (exhaustively)

% load spiking data of 15 neurons
load example15

% randomly divide it into a training set and a test set (so we can verify how well we trained)
[ncells,nsamples] = size(spikes15);
idx_train = randperm(nsamples,ceil(nsamples/2));
idx_test = setdiff(1:nsamples,idx_train);
samples_train = spikes15(:,idx_train);
samples_test = spikes15(:,idx_test);

% create a k-ising model (pairwise maxent with synchrony constraints)
model = createModel(ncells,'kising');

% train the model to a threshold of one standard deviation from the error of computing the marginals.
% because the distribution is relatively small (15 dimensions) we can explicitly represent all 2^15 states 
% in memory and train the model in an exhaustive fashion.
model = trainModel(model,samples_train,'threshold',1);

% now check the kullback-leibler divergence between the model predictions and the pattern counts in the test-set 
empirical_distribution = getEmpiricalModel(samples_test);
model_logprobs = getLogProbability(model,empirical_distribution.words);
test_dkl = dkl(empirical_distribution.logprobs,model_logprobs);
fprintf('Kullback-Leibler divergence from test set: %f\n',test_dkl);

model_entropy = getEntropy(model);
fprintf('Model entropy: %.03f   empirical dataset entropy: %.03f\n', getEntropy(model), empirical_distribution.entropy);

% get the marginals (firing rates and correlations) of the test data and see how they compare to the model predictions
marginals_data = getEmpiricalMarginals(samples_test,model);
marginals_model = getMarginals(model);

% plot them on a log scale
figure
loglog(marginals_data,marginals_model,'b*');
hold on;
minval = min([marginals_data(marginals_data>0)]);
plot([minval 1],[minval 1],'-r'); % identity line
xlabel('empirical marginal');
ylabel('predicted marginal');
title(sprintf('marginals in %d cells',ncells));


%% part 2: working with larger distributions of neurons (MCMC)

load example50

% randomly divide into train/test sets
[ncells,nsamples] = size(spikes50);
idx_train = randperm(nsamples,ceil(nsamples/2));
idx_test = setdiff(1:nsamples,idx_train);
samples_train = spikes50(:,idx_train);
samples_test = spikes50(:,idx_test);

% create an ising model (pairwise maxent)
model = createModel(50,'ising');

% train the model to a threshold of 1.5 standard deviations from the error of computing the marginals.
% because the distribution is larger (50 dimensions) we cannot explicitly iterate over all 5^20 states 
% in memory and will use markov chain monte carlo (MCMC) methods to obtain an approximation
model = trainModel(model,samples_train,'threshold',1.5);


% get the marginals (firing rates and correlations) of the test data and see how they compare to the model predictions.
% here the model marginals could not be computed exactly so they will be estimated using monte-carlo. We specify the
% number of samples we use so that their estimation will have the same amoutn noise as the empirical marginal values
marginals_data = getEmpiricalMarginals(samples_test,model);
marginals_model = getMarginals(model,'nsamples',size(samples_test,2));

% plot them on a log scale
figure
loglog(marginals_data,marginals_model,'b*');
hold on;
minval = min([marginals_data(marginals_data>0)]);
plot([minval 1],[minval 1],'-r'); % identity line
xlabel('empirical marginal');
ylabel('predicted marginal');
title(sprintf('marginals in %d cells',ncells));

% the model that the MCMC solver returns is not normalized. If we want to compare the predicted and actual probabilities
% of individual firing patterns, we will need to first normalize the model. We will use the wang-landau algorithm for
% this
disp('Normalizing model...');
model = wangLandau(model);

% the normalization factor was added to the model structure. Now that we have a normalized model, we'll use it to
% predict the frequency of activity patterns. We will start by observing all the patterns that repeated at least twice
% (because a pattern that repeated at least once may grossly overrepresent its probability and is not meaningful in this
% sort of analysis)
limited_empirical_distribution = getEmpiricalModel(samples_test,'min_count',2);


% get the model predictions for these patterns
model_logprobs = getLogProbability(model,limited_empirical_distribution.words);

% nplot on a log scale
figure
plot(limited_empirical_distribution.logprobs,model_logprobs,'bo');
hold on;
minval = min(limited_empirical_distribution.logprobs);
plot([minval 0],[minval 0],'-r');  % identity line
xlabel('empirical pattern log frequency');
ylabel('predicted pattern log frequency');
title(sprintf('activity pattern frequency in %d cells',ncells));



% Wang-landau also approximated the model entropy, let's compare it to the entropy of the empirical dataset.
% for this we want to look at the entire set, not just the set limited repeating patterns
empirical_distribution = getEmpiricalModel(samples_test);

% it will not be surprising to see that the empirical entropy is much lower than the model, this is because the
% distribution is very undersampled
fprintf('Model entropy: %.03f bits, empirical entropy (test set): %.03f bits\n',model.entropy,empirical_distribution.entropy);

% generate samples from the distribution and compute their entropy. This should give a result which is must closer to
% the entropy of the empirical distribution...
samples_simulated = generateSamples(model,numel(idx_test));
simulated_empirical_distribution = getEmpiricalModel(samples_simulated);
fprintf('Entropy of simulated data: %.03f bits\n',simulated_empirical_distribution.entropy);



%% part 2: working with MERP distributions of neurons

% load spiking data of 15 neurons
load example15

% randomly divide it into a training set and a test set (so we can verify how well we trained)
[ncells,nsamples] = size(spikes15);
idx_train = randperm(nsamples,ceil(nsamples/2));
idx_test = setdiff(1:nsamples,idx_train);
samples_train = spikes15(:,idx_train);
samples_test = spikes15(:,idx_test);

% create a MERP model with default settings
model = createModel(ncells,'merp');

% train the model to a threshold of one standard deviation from the error of computing the marginals.
% because the distribution is relatively small (15 dimensions) we can explicitly represent all 2^15 states 
% in memory and train the model in an exhaustive fashion.
model = trainModel(model,samples_train,'threshold',1);

% now check the kullback-leibler divergence between the model predictions and the pattern counts in the test-set 
empirical_distribution = getEmpiricalModel(samples_test);
model_logprobs = getLogProbability(model,empirical_distribution.words);
test_dkl = dkl(empirical_distribution.logprobs,model_logprobs);
fprintf('Kullback-Leibler divergence from test set: %f\n',test_dkl);

% create a MERP model with a specified number of projections and specified sparsity
model = createModel(ncells,'merp','nprojections',300,'sparsity',0.5);

% train the model
model = trainModel(model,samples_train,'threshold',1);

% now check the kullback-leibler divergence between the model predictions and the pattern counts in the test-set 
empirical_distribution = getEmpiricalModel(samples_test);
model_logprobs = getLogProbability(model,empirical_distribution.words);
test_dkl = dkl(empirical_distribution.logprobs,model_logprobs);
fprintf('Kullback-Leibler divergence from test set: %f\n',test_dkl);


% normalize the MERP model threshold so that the moment-generating function will give us a requested average value:
model = normalizeMerpThreshold(model,samples_train);

% retrain it
model = trainModel(model,samples_train,'threshold',1);

% now check the kullback-leibler divergence between the model predictions and the pattern counts in the test-set 
empirical_distribution = getEmpiricalModel(samples_test);
model_logprobs = getLogProbability(model,empirical_distribution.words);
test_dkl = dkl(empirical_distribution.logprobs,model_logprobs);
fprintf('Kullback-Leibler divergence from test set: %f\n',test_dkl);



