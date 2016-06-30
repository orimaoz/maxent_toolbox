% Trains a maximum-entropy model on a set of samples.
%
% Usage:
%   model_out = trainModel(model,samples)
%   model_out = trainModel(model,samples,Name,Value,...)
%
% The function will automatically choose, based on the number of dimensions in the input distribution, whether to
% compute the distribution exhaustively (in which case it will return a normalized model) or to compute the distribution
% using Monte-Carlo Markov Chain (MCMC) techniques (in which case it will return an un-normalized model, which can 
% later be approximately normalized using other functions in this package). The user can force either an exhaustive
% solution or an MCMC solution by supplying optional arugments.
% If the model is an independent model, it will extract the model parameters analytically.
%
% Arguments (mandatory):
%   model   - Maximum entropy model as returned by the createModel() function. 
%   samples - Set of samples to train the model on, in the format (ncells x nsamples). 
%             If the input dimensionality is small enough for an exhaustive computation, a Boltzmann distribution (in the
%             same format as returned by createModel) may be inputed instead of a raster, in this case the target
%             marginals will be computed in an exact manner from the distribution.   
%
% Optional arguments (in the form of Name,Value pairs):
%   threshold        - convergence threshold
%   silent           - don't print anything
%   max_steps        - limit to a maximum number of steps 
%   use_acceleration - true/false to use accelerated gradient descent (enabled by default)
%   force_exhaustive - set to true to force an exhaustive numerical solution. This entails storing in memory all
%                      2^ncells states of the distribution which can grow really fast, so is not recommended for 
%                      inputs of more than 30 cells.
%   force_mcmc -       set to true to force an MCMC solver even when a btter option is available.
%   savefile         - will constantly save the state in this file, and try to resume from it if it already exists.
%
% Output:
%   model - trained ME model. If the input size was small enough, the model will also be normalized.
%
%
% Example usage:
%   model_out = trainModel(model,raster,'threshold',0.005,'savefile','training_reentry.mat');
%
% Last update: Ori Maoz 30/06/2016
function model_out = trainModel(input_model,raster,varargin)

MAXIMUM_NCELLS_FOR_EXHAUSTIVE = 25;

if (nargin<2)
    error('Usage: model_out = trainModel(model,samples,optional_arguments)');
end

% before checking the population size, check if we got as input an empirical set or a model
if (isstruct(raster))
    % we got a model (that we can sample from)
    ncells = raster.ncells;
else
    % our input is a raster
    ncells = size(raster,1);
end

if (input_model.ncells ~= ncells)
    error('dimensionality mismatch - model and data have different number of cells');        
end

% check if the user wants to force any type of solver
p = inputParser;
addOptional(p,'force_exhaustive',false,@islogical);
addOptional(p,'force_mcmc',false,@islogical);
p.KeepUnmatched = true;
p.parse(varargin{:});


if (p.Results.force_exhaustive)
    % if the user demands an exhaustive solution, give it to him
    model_out = trainModelExhaustive(input_model,raster,varargin{:});
elseif (p.Results.force_mcmc)
    % if the user demands an MCMC solution, give it to him
    model_out = trainModelMCMC(input_model,raster,varargin{:});
elseif (strcmp(input_model.type,'indep'))
    % if this is an independent model, send it to a separate trainer because there is an analytic solution
    model_out = trainIndepModel(input_model,raster);    
else
    % check if there are few cells in the input raster so we can do an exhaustive training
    if (ncells <= MAXIMUM_NCELLS_FOR_EXHAUSTIVE)
        model_out = trainModelExhaustive(input_model,raster,varargin{:});
    else
        model_out = trainModelMCMC(input_model,raster,varargin{:});
    end
end

end