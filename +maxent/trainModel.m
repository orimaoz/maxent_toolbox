% Trains a maximum-entropy model on a set of samples.
%
% Usage:
%   [model_out, bConverged] = trainModel(model,samples)
%   [model_out, bConverged] = trainModel(model,samples,Name,Value,...)
%
% Trains a maximum entropy model on empirical data. The function will automatically choose, based on the number of 
% dimensions in the input distribution, which of two modes of operation to use:
%
% * For a small number (default <= 25) of input dimensions it will compute an exact solution and return 
%   a normalized probability distribution.
% * For a large number (default > 25) of input dimensions it will compute an approximate solution using 
%   Monte-Carlo Markov Chain (MCMC) methods and return a non-normalized distribution. This distribution 
%   can later be normalized using other functions in the toolbox such as wangLandau.
% * If the input model is an independent model, the function will return a normalized probability distribution 
%   regardless of the input dimension.
%
% The user can force either an exhaustive solution or an MCMC solution by supplying optional arguments.
%
% Arguments (mandatory):
%   model   - Maximum entropy model as returned by the createModel() function. 
%   samples - Set of samples to train the model on, in the format (ncells x nsamples). 
%             If the input dimensionality is small enough for an exhaustive computation, a Boltzmann distribution (in the
%             same format as returned by createModel) may be inputed instead of a raster, in this case the target
%             marginals will be computed in an exact manner from the distribution.   
%
% Optional arguments (in the form of Name,Value pairs):
%   threshold        - convergence threshold in units of standard deviations in the marginals. This standard deviation
%                      is estimated using the number of samples the marginal was computed from, which means that larger 
%                      datasets will be assigned tighter (actual) thresholds. 
%                      For small groups of inputs (where the model can be built in an exhaustive manner) specifying a
%                      threshold of zero will set the threshold to the quantization error of the marginals, which is
%                      equal to (1/(nsamples*2)).
%   silent           - don't print anything
%   max_steps        - limit to a maximum number of steps 
%   use_acceleration - true/false to use accelerated gradient descent (enabled by default)
%   force_exhaustive - set to true to force an exhaustive numerical solution. This entails storing in memory all
%                      2^ncells states of the distribution which can grow really fast, so is not recommended for 
%                      inputs of more than 25 cells.
%   force_mcmc       - set to true to force an MCMC solver even when a btter option is available.
%   savefile         - will constantly save the state in this file, and try to resume from it if it already exists.
%   save_delay       - delay between saves (in seconds).
%
% Output:
%   model      - trained ME model. If the input size was small enough, the model will also be normalized.
%   bConverged - true if the training process converged, false if it reached the maximum #iterations
%
%
% Example usage:
%   model = maxent.createModel(20,'pairwise');
%   model_out = maxent.trainModel(model,samples,'threshold',1,'savefile','training_reentry.mat');
%
% Last update: Ori Maoz 09/11/2017
function [model_out, bConverged] = trainModel(input_model,samples,varargin)

MAXIMUM_NCELLS_FOR_EXHAUSTIVE = 25;

if (nargin<2)
    error('Usage: model_out = trainModel(model,samples,[optional_arguments])');
end

% before checking the population size, check if we got as input an empirical set or a model
if (isstruct(samples))
    % we got a model (that we can sample from)
    ncells = samples.ncells;
else
    % our input is a raster
    ncells = size(samples,1);
end

if (input_model.ncells ~= ncells)
    error('dimensionality mismatch - model and data have different number of cells');        
end

% check if the user wants to force any type of solver
p = inputParser;
addOptional(p,'force_exhaustive',false,@islogical);   % force exhaustive solver
addOptional(p,'force_mcmc',false,@islogical);         % force MCMC solver
p.KeepUnmatched = true;
p.parse(varargin{:});


if (p.Results.force_exhaustive)
    % if the user demands an exhaustive solution, give it to him
    [model_out,bConverged] = maxent.trainModelExhaustive(input_model,samples,varargin{:});
elseif (p.Results.force_mcmc)
    % if the user demands an MCMC solution, give it to him
    [model_out,bConverged] = maxent.trainModelMCMC(input_model,samples,varargin{:});
elseif (strcmp(input_model.type,'indep'))
    % if this is an independent model, send it to a separate trainer because there is an analytic solution
    model_out = maxent.trainModelIndependent(input_model,samples);    
    bConverged = true;  % training an independent model is analytic 
else
    % check if there are few cells in the input raster so we can do an exhaustive training
    if (ncells <= MAXIMUM_NCELLS_FOR_EXHAUSTIVE)
        [model_out,bConverged] = maxent.trainModelExhaustive(input_model,samples,varargin{:});
    else
        [model_out,bConverged] = maxent.trainModelMCMC(input_model,samples,varargin{:});
    end
end

end