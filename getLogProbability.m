% Returns the log probability (in natural logarithm) of samples according to a trained maximum entropy model
%
% Usage:
%   logprobs = getLogProbability(model, samples)
%   logprobs = getLogProbability(model, samples,Name,Value,...)
%
%
% The function accepts a model and an array of samples and returns the log probability for each sample 
% according to the model. If the model has not been normalize, the results will also be un-normalized.
%
% Arguments (mandatory):
%   model   - Maximum entropy model as returned by trainModel()
%   samples - Set of samples to train the model on, in the format (ncells x nsamples). 
%
% Optional arguments (in the form of Name,Value pairs):
%   normalize   - force with or without normalization. If this value is true but the model is un-normalized, the
%                 function will issue a warning.
%
% Output:
%   logprobs - an (1xnsamples) vector of log probabilities in natural logarithm.
%
%
% Last update: Ori Maoz 30/06/2016
function logprobs = getLogProbability(model, samples,varargin)

p = inputParser;
addOptional(p,'normalize',nan);    % silent mode - don't print anything

p.parse(varargin{:});
normalize = p.Results.normalize;

% if the user specified something as to normalization
if ~isnan(normalize)
    
    % the user specifically asked for normalization
    if (normalize)
        if (~isfield(model,'z'))
            warning('Model is un-normalized!');
        end        
    else    
        % the user specifically asked for no normalization
        if (isfield(model,'z'))
            model = rmfield(model,'z');
        end
    end
end
    
% Now delegate to the MEX implementation
logprobs = mexLogProbability(samples,model);

end