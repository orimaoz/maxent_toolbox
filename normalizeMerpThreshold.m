% Tunes the threholds of a MERP projection so that its marginals will be centered around some value.
%
% Usage:
%   model_out = normalizeMerpThreshold(model,samples)
%   model_out = normalizeMerpThreshold(model,samples,Name,Value,...)
%
% This function will change the thresholds used by a MERP model so that the average activity of the hidden layer
% will be around the requested number. It can either set all the thresholds to the same value (in which case it will
% try to match the mean activity to the requested value) or each of the values individually.
%
% Arguments (mandatory):
%   model   - MERP model as returned by the createModel() function. 
%   samples - Set of samples used to normalize the model, in the format (ncells x nsamples). 
%
% Optional arguments (in the form of Name,Value pairs):
%   target           - target marginal value (default 0.3)
%   individual       - set to true if you want to tune each threhold by itself, false to make all the threholds the
%                      same value and only look at the average marginal (default false)
%   threshold_range  - vector of [min_thresh, max_thresh], smallest and largest possible thresholds
%   silent           - don't print anything
%
% Output:
%   model - trained ME model. If the input size was small enough, the model will also be normalized.
%
%
% Example usage:
%   model_out = normalizeMerpThreshold(model,samples,'target',0.2)
%
% Last update: Ori Maoz 03/07/2016
%
function model_out = normalizeMerpThreshold(model,samples,varargin)

DEFAULT_TARGET = 0.3;
DEFAULT_INDIVIDUAL = false;
DEFAULT_SILENT = false;

NUMBER_OF_STEPS = 20;

% this is the value range in which we will tune it 
DEFAULT_THRESHOLD_MAX = 50;
DEFAULT_THRESHOLD_MIN = 0.000001;

if ~(strcmpi(model.type,'merp') || strcmpi(model.type,'merpsparse'))
    error('Model must be a MERP model');
end


% set up input arguments;
p = inputParser;
addOptional(p,'target',DEFAULT_TARGET,@isnumeric); % default target rate
addOptional(p,'individual',DEFAULT_INDIVIDUAL,@islogical);    % true to tune each projection separately
addOptional(p,'threshold_range',[DEFAULT_THRESHOLD_MIN,DEFAULT_THRESHOLD_MAX]); % range of possible thresholds
addOptional(p,'silent',DEFAULT_SILENT,@islogical); % don't print anything
p.parse(varargin{:});
target_rate = p.Results.target;
tune_individual = p.Results.individual;
threshold_min = p.Results.threshold_range(1);
threshold_max = p.Results.threshold_range(2);


nfactors = numel(model.factors);

% initial thresholds - start from the middle of the inputed range and converge from there
upper_values = ones(size(model.factors)) * threshold_max;
lower_values = ones(size(model.factors)) * threshold_min;
model.threshold = (upper_values + lower_values)/2;

% use binary search to zoom in on the new values
for i = 1:NUMBER_OF_STEPS
       
    % find the current rates
    rates = mexGetMarginals(samples,model);
    
       
    if (tune_individual)
        
        % set each change direction separately
        idx_too_low = rates > target_rate;
        idx_too_high = rates < target_rate;
    else
        % make them all change in the same direction
        if (median(rates) > target_rate)
            idx_too_low = 1:nfactors;
            idx_too_high = [];            
        else
            idx_too_high = 1:nfactors;
            idx_too_low = [];            
        end
    end
        
    % now limit our range and modify the model
    upper_values(idx_too_high) = (upper_values(idx_too_high) + lower_values(idx_too_high))/2;
    lower_values(idx_too_low) = (upper_values(idx_too_low) + lower_values(idx_too_low))/2;
    model.threshold = (upper_values + lower_values)/2;
    
    fprintf('.')
         
end

model_out = model;




end

