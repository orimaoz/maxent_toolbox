% Trains a maximum-entropy model on a set of samples using Monte Carlo Markov Chain (MCMC). This is a sub-function 
% invoked by trainModel() and typically should not be invoked directly.
%
% Usage:
%   model_out = trainModelMCMC(model,samples)
%   model_out = trainModelMCMC(model,samples,Name,Value,...)
%
% Input:
%   refer to trainModel.m
%
% Output:
%   model_out - trained ME model, un-normalized.
%
% Last update: Ori Maoz 30/06/2016
function [model] = trainModelMCMC(input_model,raster,varargin)


last_save = tic;

% how many recent iterations we use for the distribution of marginal errors
GRADIENT_MEMORY_SIZE = 50;

% how many initial steps of GIS to use before switching back to gradient descent
GIS_INITIAL_STEPS = 50;

% how many samples we sample in the beginning, and how quickly we raise the amount 
DEFAULT_MIN_NSAMPLES_BASE = 100;
DEFAULT_NSAMPLES_INCREASE = 1.01;

% how many previous gradients we store, we use this to estimate the changes in step size
GRADIENT_NORM_COUNT = 10;


% number of seconds after we save our current status to a file
DEFAULT_TIME_BETWEEN_SAVES = 60;


DEFAULT_NUM_STEPS = inf;
DEFAULT_MAX_ERROR_THRESHOLD = 1.3;


% check dimensionality
[ncells, data_nsamples] = size(raster);
if (data_nsamples < ncells)    
    warning('Input raster must be of the form (nsamples x ncells), are you sure you did not transpose it?');
end


% parse our optional arguments
p = inputParser;
addOptional(p,'threshold',DEFAULT_MAX_ERROR_THRESHOLD,@isnumeric); % convergence threshold
addOptional(p,'savefile','');                 % save file for re-entrant code
addOptional(p,'max_steps',DEFAULT_NUM_STEPS); % maximum number of steps that we run
addOptional(p,'use_acceleration',true,@islogical);    % true if we want to use accelerated gradient descent (nesterov)
addOptional(p,'save_delay',DEFAULT_TIME_BETWEEN_SAVES);    % time between re-entry file saves
addOptional(p,'max_nsamples',nan,@isnumeric);    % maximum number of samples used to estimate marginals
addOptional(p,'nsamples_increase',DEFAULT_NSAMPLES_INCREASE,@isnumeric);    %  how quicly we raise the number of samples
addOptional(p,'silent',false,@islogical);    % silent mode - don't print anything


p.parse(varargin{:});
reentry_filename = p.Results.savefile;
num_steps = p.Results.max_steps;
use_acceleration = p.Results.use_acceleration;
time_between_saves = p.Results.save_delay;
max_nsamples = p.Results.max_nsamples;
requested_threshold = p.Results.threshold;
nsamples_increase = p.Results.nsamples_increase;
silent = p.Results.silent;


if isnan(max_nsamples)
    % if the user did not specify how many samples to use in the MCMC simulation, choose a number which is large
    % enough to obtain a reliable measurement for the requested confidence interval
    max_nsamples = ceil(2*data_nsamples /((requested_threshold)^2));
end

% The internal threshold is a function of both the variance of the MCMC (reflected in max_nsamples) and the requeted
% threshold. Choose an internal threshold that would force the distance to the target to be within the specified range.
max_error_threshold = sqrt(requested_threshold^2 + data_nsamples/max_nsamples);

internal_print(sprintf('Training to threshold: %.03f standard deviations',requested_threshold));
internal_print('Maximum samples: %d   maximum MSE: %.03f',max_nsamples,max_error_threshold^2);

% make sure that the reentry filename is well-structured
if ~isempty(reentry_filename)
    [pathstr,name,~] = fileparts(reentry_filename);
    reentry_filename  = fullfile(pathstr,[name '.mat']);
    
    internal_print('Using save file %s for re-entry\n',reentry_filename);
end

stderr=nan(1,numel(input_model.factors));  % square errors for each factor


%if the re-entry file exists, load running status from it
i=0;
if exist(reentry_filename, 'file') == 2
    load(reentry_filename);
    
    % reset the time of the last save so we won't re-save immediately
    last_save = tic;

else
    
    % if it does not exist, initialize all that should be initialized
    model = input_model;


    % the diameter parameter multipliers the step size so that the convergence could be adapted to models of different
    % types. If it seems that the solver cannot converge because the step size is too large, decreasing model.D might help.
    if isfield(model,'step_scaling')
        if (numel(model.step_scaling) ~= numel(model.factors))
            error('model step scaling does not match number of factors');
        end
        step_scaling = model.step_scaling;
    else
        warning('No model step scaling provided, using default');
        step_scaling = ones(size(model.factors));
    end


    % get the empirical values of the marginals - our goal is to fit the model marginals to these.
    empirical_marginals = mexEmpiricalMarginals(raster,model);


    % estimate the standard deviation of each of the marginals to see how close we should get to it.
    [phat marginals_std] = binofit(round(empirical_marginals * data_nsamples),data_nsamples,0.32);
    empirical_marginals_std = max(abs(repmat(phat',[1 2]) - marginals_std),[],2)';    
    
    % number of largrange multipliers of the model
    nfactors = numel(model.factors);

    % here we will keep the recent gradient norm history, we will use this to estimate the step size
    gradient_norm_history = ones(1,GRADIENT_NORM_COUNT);
    gradient_norm_idx = 1;

    % here we will keep the absolute deviations of marginals, we will use this to check for convergence
    gradient_memory = nan(GRADIENT_MEMORY_SIZE,nfactors);
    gradient_memory_idx = 1;

    prev_nsamples = DEFAULT_MIN_NSAMPLES_BASE;

    if (use_acceleration)
        nesterov_acceleration = Nesterov;  % initialize Nesterov accelerated gradient descent
        y = model.factors;
    end

    % Starting point of the markov chain. This variable will be changed by the gibbs sampler so that each markov
    % chain will continue from the end of the previous chain.
    x0 = uint32(zeros(1,ncells));
end

last_print_time = tic;
while i < num_steps
          
       
    i = i + 1;



    % limit the growth speed of the # of samples
    nsamples = ceil(prev_nsamples * nsamples_increase);
    nsamples = min(nsamples,max_nsamples);
    
 
    % generate samples and compute marginals from them
    params.separation = 1;
    params.burnin = 0;
    x = (gibbsSampler(x0,nsamples,model,params));
    model_marginals = mexEmpiricalMarginals(x,model);         
    
    prev_nsamples = nsamples;
    
    % factor_gradient of the entropy function
    factor_gradient = model_marginals - empirical_marginals;
        
    % save the gradient, which is also how far we are from the target marginals
    gradient_memory(gradient_memory_idx,:) = abs(factor_gradient);    
    gradient_memory_idx = mod(gradient_memory_idx,GRADIENT_MEMORY_SIZE)+1;    
   
   
    % compute the variance of estimation error for each marginal across time, normalized by the expected standard
    % deviation of this marginal
    if (i >= GRADIENT_MEMORY_SIZE)
        stderr = sqrt(mean(gradient_memory.^2)) ./ empirical_marginals_std;
    end    
    mean_error = mean(stderr);
    [max_error, max_idx] = max(stderr);
    
    
    % Show difference of firing rates and DKL to target (only print once every maximum 1 second)        
    t=toc(last_print_time);
    if t > 1
        %internal_print('%02d/%02d samples=%d MSE mean:%.04f max: %.04f [%d]',i,num_steps,nsamples,mean_error^2,max_error^2,max_idx);
        internal_print('%02d/%02d samples=%d  MSE=%.03f (mean), %.03f (max) [%d]',i,num_steps,nsamples,mean_error^2,max_error^2,max_idx);
        last_print_time = tic;
    end


    if (max_error < max_error_threshold)
        internal_print('converged (marginals match)');
        break;
    end

    
    % remember the norm of the gradient.
    % an estimation of this will be used to compute the step size
    gradient_norm_history(gradient_norm_idx) = dot(factor_gradient,factor_gradient);
    gradient_norm_idx = mod(gradient_norm_idx,GRADIENT_NORM_COUNT)+1;
    
    
    % Step size computation:
    % from theoretical SGD analysis, we know that if our gradient is noisy and:
    % E[||gradient||^2] < G^2
    % then we get good regret guarantees when we choose
    % step_size = D / (G * sqrt(t))    (t = current step index)
    % D is an upper bound on the diameter of the problem and is therefore problem-specific. 
    G2 = mean(gradient_norm_history);
    curr_step_size = step_scaling ./ sqrt(G2 * i);    

    
    % for a certain number of initial steps we will perform Generalized Iterative Scaling which is more stable than
    % gradient descent when we are still far from convergence
    if i < GIS_INITIAL_STEPS        
        
        %GIS update
        C = min(1.1*sum(empirical_marginals),nfactors);        
        delta = (log(empirical_marginals ./ model_marginals)) / C;
        delta(~isfinite(delta))=0;
        model.factors = model.factors - delta;            
        
        % save momentum for accelerated gradient descent (for when we switch)
        y = model.factors + curr_step_size .* factor_gradient;
    else        
        % use gradient descent        
        if (use_acceleration)
            
            % Use Nesterov's accelerated gradient ascent
            nesterov_acceleration = nesterov_acceleration.nextGamma;
            gamma = nesterov_acceleration.gamma;
            
            yprev = y;
            y = model.factors + curr_step_size .* factor_gradient;
            model.factors = (1-gamma) .* y + gamma .* yprev;
        else       
            
            % Vanilla gradient descent
            model.factors = model.factors + curr_step_size .* factor_gradient;
        end
    end
    
    % check if we should save re-entry data
    if ~isempty(reentry_filename)
        time_from_last_save = toc(last_save);
        if (time_from_last_save > time_between_saves)
            internal_print('saving re-entry file...');
            save(reentry_filename);
            last_save = tic;
        end
    end
    
end

if (i==num_steps)
    internal_print('Reached maximum iterations, stopping.');
end


    % print a message only if message printing has not been disabled ("silent mode")
    function internal_print(varargin)
        if (~silent)
            disp(sprintf(varargin{:}));
        end        
    end
    
end

