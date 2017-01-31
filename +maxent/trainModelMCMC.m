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
%   bConverged - true if the training process converged, false if it reached the maximum #iterations
%
% Last update: Ori Maoz 25/01/2017
function [model, bConverged] = trainModelMCMC(input_model,raster,varargin)


last_save = tic;

% how many recent iterations we use for the distribution of marginal errors
GRADIENT_MEMORY_SIZE_INCREMENT = 50;
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


% number of steps after which we start by increase the sample size
NSTEPS_ADJUST_SAMPLE_SIZE = 5000;

% MSE for which we decide that we have diverged and should reset the training
DIVERGENCE_THRESHOLD = 100;


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
original_max_samples = max_nsamples;


persistent_state.model = input_model;

% number of largrange multipliers of the model
nfactors = maxent.getNumFactors(persistent_state.model);
   


% make sure that the reentry filename is well-structured
if ~isempty(reentry_filename)
    [pathstr,name,~] = fileparts(reentry_filename);
    reentry_filename  = fullfile(pathstr,[name '.mat']);
    
    internal_print('Using save file %s for re-entry\n',reentry_filename);
end



%if the re-entry file exists, load running status from it
if exist(reentry_filename, 'file') == 2
    load(reentry_filename);
    
    % reset the time of the last save so we won't re-save immediately
    last_save = tic;

else
    
    persistent_state.i=0;
    persistent_state.curr_iteration = 0;
    persistent_state.time_since_last_adjustment = 0;
    
    persistent_state.max_nsamples = max_nsamples;
    
    % if it does not exist, initialize all that should be initialized
    persistent_state.model = input_model;
    
    persistent_state.stderr=nan(1,nfactors);  % square errors for each factor
    
    % here we will keep the recent gradient norm history, we will use this to estimate the step size
    persistent_state.gradient_norm_history = ones(1,GRADIENT_NORM_COUNT);
    persistent_state.gradient_norm_idx = 1;

    % here we will keep the absolute deviations of marginals, we will use this to check for convergence
    persistent_state.gradient_memory = nan(GRADIENT_MEMORY_SIZE,nfactors);
    persistent_state.gradient_memory_idx = 1;
    
    persistent_state.prev_nsamples = DEFAULT_MIN_NSAMPLES_BASE;

    if (use_acceleration)
        persistent_state.nesterov_acceleration = maxent.Nesterov;  % initialize Nesterov accelerated gradient descent
        persistent_state.y = maxent.getFactors(persistent_state.model);
    end

    % Starting point of the markov chain. This variable will be changed by the gibbs sampler so that each markov
    % chain will continue from the end of the previous chain.
    persistent_state.x0 = uint32(zeros(1,ncells));

    
end
    



% the diameter parameter multipliers the step size so that the convergence could be adapted to models of different
% types. If it seems that the solver cannot converge because the step size is too large, decreasing model.D might help.

step_scaling = maxent.getStepScaling(persistent_state.model);
if ~isnan(step_scaling)
    if (numel(step_scaling) ~= nfactors)
        error('model step scaling does not match number of factors');
    end
else
    warning('No model step scaling provided, using default');
    step_scaling = ones(1,nfactors);
end



% get the empirical values of the marginals - our goal is to fit the model marginals to these.
empirical_marginals = maxent.getEmpiricalMarginals(raster,persistent_state.model);


% estimate the standard deviation of each of the marginals to see how close we should get to it.
[phat marginals_std] = binofit(round(empirical_marginals * data_nsamples),data_nsamples,0.32);
empirical_marginals_std = max(abs(repmat(phat',[1 2]) - marginals_std),[],2)';    
    

% The internal threshold is a function of both the variance of the MCMC (reflected in max_nsamples) and the requeted
% threshold. Choose an internal threshold that would force the distance to the target to be within the specified range.
max_error_threshold = sqrt(requested_threshold^2 + data_nsamples/persistent_state.max_nsamples);


internal_print(sprintf('Training to threshold: %.03f standard deviations',requested_threshold));
internal_print('Maximum samples: %d   maximum MSE: %.03f',max_nsamples,max_error_threshold^2);

last_print_time = tic;
while persistent_state.i < num_steps
          
       
    persistent_state.i = persistent_state.i + 1;  % steps in the main loop
    persistent_state.curr_iteration = persistent_state.curr_iteration + 1;   % steps for the solver


    % if we seem to be running a very long time, increase the sample size
    persistent_state.time_since_last_adjustment = persistent_state.time_since_last_adjustment + 1;
    if persistent_state.time_since_last_adjustment > NSTEPS_ADJUST_SAMPLE_SIZE
        persistent_state.time_since_last_adjustment = 0;
        
        
        % every time increase the number of samples by one "jump"
        persistent_state.max_nsamples = persistent_state.max_nsamples + original_max_samples;

        % also remember a longer history of steps for a more accurate variance estimation
        GRADIENT_MEMORY_SIZE = GRADIENT_MEMORY_SIZE + GRADIENT_MEMORY_SIZE_INCREMENT;
        
        % adjust threshold accordingly to account for the more accurate measurement
        max_error_threshold = sqrt(requested_threshold^2 + data_nsamples/persistent_state.max_nsamples);

        % apprently we are way off, reset back to initial guess and let it reconverge with stricter sampling
        if (mean_error^2 > DIVERGENCE_THRESHOLD)
            internal_print('Solver diverged, resetting parameters and increasing max sample size (maximum MSE: %.03f)',max_error_threshold); 

            
            % reset everything- factors, history etc
            
            model = maxent.setFactors(model,zeros(1,maxent.getNumFactors(model)));
            
            % here we will keep the recent gradient norm history, we will use this to estimate the step size.
            % we are setting it such that it overestimates because we want to be on the careful side
            % and start out with smaller steps. It will gradually be overridden with real data as we proceed.
            persistent_state.gradient_norm_history = ones(1,GRADIENT_NORM_COUNT) * nfactors;
            persistent_state.gradient_norm_idx = 1;

            % here we will keep the absolute deviations of marginals, we will use this to check for convergence
            persistent_state.gradient_memory = nan(GRADIENT_MEMORY_SIZE,nfactors);
            persistent_state.gradient_memory_idx = 1;
            
%             mean_error = nan;
%             max_error = nan;
%             stderr = nan;
%             max_idx= nan;
%             
            persistent_state.curr_iteration = 0;
            
            if (use_acceleration)
                persistent_state.nesterov_acceleration = maxent.Nesterov;  % initialize Nesterov accelerated gradient descent
                persistent_state.y = maxent.getFactors(persistent_state.model);
            end
        else            
            internal_print('Unable to converge, increasing max sample size (maximum MSE: %.03f)',max_error_threshold);      
        end
        
    end
    
    

    % limit the growth speed of the # of samples
    nsamples = ceil(persistent_state.prev_nsamples * nsamples_increase);
    nsamples = min(nsamples,persistent_state.max_nsamples);
    
 
    % generate samples and compute marginals from them
    params.separation = ncells;
    params.burnin = 0;
    x = maxent.generateSamples(persistent_state.model, nsamples, 'x0',persistent_state.x0,'burnin',0,'separation',ncells);
    model_marginals = maxent.getEmpiricalMarginals(x,persistent_state.model);
    
    
    persistent_state.prev_nsamples = nsamples;
    
    % factor_gradient of the entropy function
    factor_gradient = model_marginals - empirical_marginals;
        
    % save the gradient, which is also how far we are from the target marginals
    persistent_state.gradient_memory(persistent_state.gradient_memory_idx,:) = abs(factor_gradient);    
    persistent_state.gradient_memory_idx = mod(persistent_state.gradient_memory_idx,GRADIENT_MEMORY_SIZE)+1;    
   
   
    % compute the variance of estimation error for each marginal across time, normalized by the expected standard
    % deviation of this marginal
    if (persistent_state.curr_iteration >= GRADIENT_MEMORY_SIZE)
        stderr = sqrt(mean(persistent_state.gradient_memory.^2)) ./ empirical_marginals_std;
    else
        stderr = nan;
    end    
    mean_error = mean(stderr);
    [max_error, max_idx] = max(stderr);
    
    
    % Show difference of firing rates and DKL to target (only print once every maximum 1 second)        
    t=toc(last_print_time);
    if t > 1
        %internal_print('%02d/%02d samples=%d MSE mean:%.04f max: %.04f [%d]',i,num_steps,nsamples,mean_error^2,max_error^2,max_idx);
        internal_print('%02d/%02d samples=%d  MSE=%.03f (mean), %.03f (max) [%d]',persistent_state.i,num_steps,nsamples,mean_error^2,max_error^2,max_idx);
        last_print_time = tic;
    end


    if (max_error < max_error_threshold)
        internal_print('converged (marginals match)');
        bConverged = true;
        break;
    end

    
    % remember the norm of the gradient.
    % an estimation of this will be used to compute the step size
    persistent_state.gradient_norm_history(persistent_state.gradient_norm_idx) = dot(factor_gradient,factor_gradient);
    persistent_state.gradient_norm_idx = mod(persistent_state.gradient_norm_idx,GRADIENT_NORM_COUNT)+1;
    
    
    % Step size computation:
    % from theoretical SGD analysis, we know that if our gradient is noisy and:
    % E[||gradient||^2] < G^2
    % then we get good regret guarantees when we choose
    % step_size = D / (G * sqrt(t))    (t = current step index)
    % D is an upper bound on the diameter of the problem and is therefore problem-specific. 
    G2 = mean(persistent_state.gradient_norm_history);
    curr_step_size = step_scaling ./ sqrt(G2 * persistent_state.curr_iteration);    

    
    % for a certain number of initial steps we will perform Generalized Iterative Scaling which is more stable than
    % gradient descent when we are still far from convergence
    if persistent_state.curr_iteration < GIS_INITIAL_STEPS        
        
        %GIS update
        C = min(1.1*sum(empirical_marginals),nfactors);        
        delta = (log(empirical_marginals ./ model_marginals)) / C;
        delta(~isfinite(delta))=0;
        persistent_state.model = maxent.setFactors(persistent_state.model,maxent.getFactors(persistent_state.model) - delta);
        
        % save momentum for accelerated gradient descent (for when we switch)
        persistent_state.y = maxent.getFactors(persistent_state.model)+ curr_step_size .* factor_gradient;
    else        
        % use gradient descent        
        if (use_acceleration)
            
            % Use Nesterov's accelerated gradient ascent
            persistent_state.nesterov_acceleration = persistent_state.nesterov_acceleration.nextGamma;
            gamma = persistent_state.nesterov_acceleration.gamma;
            
            persistent_state.yprev = persistent_state.y;
            persistent_state.y = maxent.getFactors(persistent_state.model)+ curr_step_size .* factor_gradient;
            persistent_state.model = maxent.setFactors(persistent_state.model,(1-gamma) .* persistent_state.y + gamma .* persistent_state.yprev);
        else       
            
            % Vanilla gradient descent
            persistent_state.model = maxent.setFactors(persistent_state.model,maxent.getFactors(persistent_state.model) + curr_step_size .* factor_gradient);

        end
    end
    
    % check if we should save re-entry data
    if ~isempty(reentry_filename)
        time_from_last_save = toc(last_save);
        if (time_from_last_save > time_between_saves)
            internal_print('saving re-entry file...');
            
            %save(reentry_filename);
            
            % save the variables necessary to resume operation if we are killed and re-started with the same input.
            save(reentry_filename,'persistent_state');
            
            
            
            last_save = tic;
        end
    end
    
end

% check if we stopped because we converged or finished the #iterations
if (persistent_state.i==num_steps)
    internal_print('Reached maximum iterations, stopping.');
    bConverged = false;
end

 
% return the solved model
model = persistent_state.model;


    % print a message only if message printing has not been disabled ("silent mode")
    function internal_print(varargin)
        if (~silent)
            disp(sprintf(varargin{:}));
        end        
    end
    
end

