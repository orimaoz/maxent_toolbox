% Trains a maximum-entropy model on a set of samples backtracking and adaptive step choice. This is a sub-function 
% invoked by trainModel() and typically should not be invoked directly.
%
% Usage:
%   model_out = trainModelExhaustive(model,samples)
%   model_out = trainModelExhaustive(model,samples,Name,Value,...)
%
% Input:
%   refer to trainModel.m
%
% Output:
%   model_out - trained ME model (normalized)
%
% Last update: Ori Maoz 30/06/2016

function model = trainModelExhaustive(input_model,raster,varargin)

if nargin<2
    error('Usage: trainModelComplete(input_model,raster,threshold)');
end

MAXIMUM_EXHAUSTIVE_NCELLS = 30;

% number of seconds after we save our current status to a file. Currently not active here.
DEFAULT_TIME_BETWEEN_SAVES = 60;


% factor_gradient ascent maximum number of steps
DEFAULT_NUM_STEPS = Inf;
CONVERGED_STEP_SIZE = 10^-8;

beta = 0.5;    % parameter for backtracking line search
alpha = 0.5;   % parameter for backtracking line search
step_size = 0.01; % initial step size for factor_gradient descent (line search modifies this)


if (isstruct(raster))
    % we got a model (that we can compute the marginals directly on)
    model_base = raster;
    ncells = model_base.ncells;
    use_exact_marginals = true;
    
else

    % we got a raster and we will estimate the marginals using it
    [ncells, nsamples] = size(raster);
    if (nsamples < ncells)

        error('Input raster must be of the form (ncells x nsamples)');
    end
    use_exact_marginals = false;

end


if (ncells > MAXIMUM_EXHAUSTIVE_NCELLS)
    error(sprintf('I refuse to compute an exhaustive model for more than %d cells. Seriously, you don''t want to try.',MAXIMUM_EXHAUSTIVE_NCELLS));
    
end

DEFAULT_L1_THRESHOLD = 1;


% parse our optional arguments
p = inputParser;
addOptional(p,'threshold',DEFAULT_L1_THRESHOLD);
addOptional(p,'savefile','');
addOptional(p,'force_complete',true);
addOptional(p,'max_steps',DEFAULT_NUM_STEPS);  % maximum number of steps that we run
addOptional(p,'observable_accuracy',1000000);  % how strict we are when matching observables
addOptional(p,'use_acceleration',true,@islogical);    % acceleration
addOptional(p,'save_delay',DEFAULT_TIME_BETWEEN_SAVES); % time between re-entry file saves. Currently not active here.
addOptional(p,'silent',false,@islogical);    % silent mode - don't print anything


p.parse(varargin{:});
threshold = p.Results.threshold;
reentry_filename = p.Results.savefile;
num_steps = p.Results.max_steps;
use_acceleration = p.Results.use_acceleration;
silent = p.Results.silent;
num_steps = p.Results.max_steps;
observable_accuracy = p.Results.observable_accuracy;



model = input_model;

if (use_exact_marginals)
    % we have the input in its analytical form so we can get an exact value for the observables
    % we will compute it later after we have the individual probabilities of all patterns
else            

end



nfactors = numel(model.factors);

backtrack_now = false;


current_dkl = inf;
previous_dkl = inf;


% set of all possible firing patterns
npatterns = (2^ncells);

if (use_exact_marginals)
    % if we have the exact observables then we need all the patterns for DKL
    unique_words = logical(de2bi(0:(npatterns-1)))';
    
    % get probabilities of the input patterns
    lp = mexLogProbability(unique_words,model_base);
    lp = lp - logsumexp(lp(:));
    empirical_probs = exp(lp);
    
    % get the observables - our goal is to fit the model marginals to these.
    empirical_marginals = mexGetMarginals(unique_words,model,empirical_probs);
    
    % compute confidence intervals for the observables. These actually are all zero (we compute them exactly) 
    % but we will still need to make some arbitrary choice when to stop. This computation makes sure that whatever
    % threshold we choose will be inforced in a stricter fashion on the marginals which are closer to the 
    % zero because the relative damage there is higher
    nsamples = observable_accuracy;
    [phat, marginals_pci] = binofit(round(empirical_marginals * nsamples),nsamples,0.32);
    empirical_marginals_std = max(abs(repmat(phat',[1 2]) - marginals_pci),[],2)';    
    
else

    % we only need the actually observed patterns in order to compute the DKL during the descent, so find out what those are
    % (to save time later).

    [unique_words I J] = unique(raster','rows','sorted');
    counts = hist(J,numel(I));  % appearance count of each word
    unique_words = uint32(unique_words');

    % normalize to probabilities
    empirical_probs = counts / sum(counts);
    
    % get the empirical values of the marginals - our goal is to fit the model marginals to these.
    empirical_marginals = mexGetMarginals(raster,model);

    % estimate the standard deviation of each of the marginals to see how close we should get to it.
    [phat, marginals_pci] = binofit(round(empirical_marginals * nsamples),nsamples,0.32);
    empirical_marginals_std = max(abs(repmat(phat',[1 2]) - marginals_pci),[],2)';    
    confidence_interval_lower = marginals_pci(:,1)';
    confidence_interval_upper = marginals_pci(:,2)';
    
end

previous_model = model;

factor_gradient = zeros(1,nfactors);

            
if (use_acceleration)
    nesterov_acceleration = Nesterov;  % initialize Nesterov accelerated gradient descent
    y = model.factors;
end

    
last_print_time = tic;
i=0;
while i < num_steps
        
    
        % get model marginals in advance, so that we will have the Z (partition function)
        % which we will require in order to compute the KL-divergence        
        [temp_model_marginals z] = mexExhaustiveMarginals(model);
        model.z = -log(z);    
        
        % compute probabilities for all possible patterns and then compute average marginals
        model_logprobs = mexLogProbability(unique_words,model);
        model_probs = exp(model_logprobs);

       
        previous_dkl = current_dkl;
        current_dkl = dkl(empirical_probs,model_probs);
        
        if isinf(current_dkl)
            %error('Received infinity DKL result');
            warning('Received infinity DKL result');
            backtrack_now = true;
        
        
        else
        
            % if our improvement was not enough, backtrack to previous step and decrease step size
            if ~(current_dkl < previous_dkl - alpha * step_size * norm(factor_gradient)^2)

                % if we are using nesterov's acceleration then we are not guaranteed to actually descend every time,
                % repeat the check with a regular gradient descent step and backtrack only if it fails too
                if (use_acceleration)


                    test_model = previous_model;
                    test_model.factors = (test_model.factors + factor_gradient * step_size);        


                    [tm z] = mexExhaustiveMarginals(test_model);
                    test_model.z = -log(z);    

                    model_logprobs = mexLogProbability(unique_words,test_model);
                    model_probs = exp(model_logprobs);
                    test_dkl = dkl(empirical_probs,model_probs);

                    % if it still fails this new test then we should really backtrack
                    if ~(test_dkl < previous_dkl - alpha * step_size * norm(factor_gradient)^2)

                        backtrack_now = true;

                    else
                        %disp('failed nesterov backtrack but passed regular check.');
                        backtrack_now = false;                                      
                    end
                else
                    % we are not using acceleration so we actually should backtrack
                    backtrack_now = true;                
                end




                if (step_size < CONVERGED_STEP_SIZE)
                    internal_print('converged.');
                    break;
                end

            else
                backtrack_now = false;            
            end
        end
            
        % if we passed the backtracking phase, get marginals and compute the next step
        if (~backtrack_now)
            i = i + 1;
            
            % it is OK to proceed, compute a new gradient for the next step            
            model_marginals = temp_model_marginals;
            
            % factor_gradient of the entropy function
            factor_gradient = model_marginals - empirical_marginals;
        else
            % back to previous step
            model = previous_model;
            current_dkl = previous_dkl;

            % decrease step size            
            step_size = step_size * beta;
            internal_print(sprintf('step size: %f',step_size));
        end                        
    
    
            
    if (step_size < CONVERGED_STEP_SIZE)
        break;
    end

    
    % estimate the error
    if (use_exact_marginals)
        marginal_errors = abs(factor_gradient);
        normalized_empirical_errors = marginal_errors;
    else
        marginal_errors = abs(factor_gradient);
        normalized_empirical_errors = marginal_errors ./ empirical_marginals_std;
    end

    if (~backtrack_now)
        % Show difference of firing rates and DKL to target
        
        [maxerror, maxidx] = max(normalized_empirical_errors);
        meanerror = mean(normalized_empirical_errors);
        t=toc(last_print_time);
        if t > 1
            internal_print(sprintf('%02d/%02d   std mean:%.04f max: %.04f [%d]   DKL: %.03f',i,num_steps,meanerror,maxerror,maxidx,current_dkl));
            last_print_time = tic;
        end
    end
    
    if (threshold > 0 && all(normalized_empirical_errors < threshold))
        % a threshold of a certain number of standard deviations (from expected mean STD) was given
        internal_print('converged (marginals match)');
        break;        
    elseif (threshold ==0 && all(marginal_errors < 1/(2*nsamples)))  
        % provided threshold was 0, we set the threshold at the quantization error
        internal_print('converged (marginals match)');
        break;        
    end

    % perform a step
    previous_model = model;
    if (use_acceleration)
                
        % Nesterov's accelerated factor_gradient ascent
        if (~backtrack_now)
            yprev = y;
            nesterov_acceleration = nesterov_acceleration.nextGamma;
        end
        
        gamma = nesterov_acceleration.gamma;
        
        y = (model.factors) + step_size * factor_gradient;
        model.factors = (1-gamma) .* y + gamma .* yprev;
                

    else
        % regular slow factor_gradient ascent
        model.factors = (model.factors + factor_gradient * step_size);
    end

end

if (i==num_steps)
    internal_print('Reached maximum iterations, stopping.');
end

internal_print(sprintf('std mean:%.04f max: %.04f [%d]   DKL: %.03f',meanerror,maxerror,maxidx,current_dkl));

% print a message only if message printing has not been disabled ("silent mode")
function internal_print(varargin)
    if (~silent)
        disp(sprintf(varargin{:}));
    end        
end


end

