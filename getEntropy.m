% Returns the empirical entropy for a dataset or computes the entropy for a model.
% This function only computes the entropy for small models (because it does so exhaustively), in order to obtain
% an approximation of the entropy of large distributions use the wangLandau function.
%
% Usage:
%   entropy = getEntropy(model)
%   entropy = getEntropy(samples)
%
% Arguments (mandatory):
%   model   - Maximum entropy model as returned by trainModel()
%   samples - an (ncells x nsamples) set of samples
%
% Output:
%   entropy - entropy in bits
%
% Last update: Ori Maoz 05/07/2016
function entropy = getEntropy(input)

MAX_EXHAUSTIVE_SIZE = 25;

% check if the input is a model or a dataset
if isfield(input,'type')
    model = input;
    if (model.ncells > MAX_EXHAUSTIVE_SIZE)
        error('model is too big to exhaustively compute the entropy, use wangLandau.m instead');
    end
    
    npatterns = 2^model.ncells;
    
    % get the probabilities of all the possible patterns
    unique_words = logical(de2bi(0:(npatterns-1)))';
    
    % get probabilities of the input patterns
    logprobs = getLogProbability(model,unique_words);
    
    % compute entropy (it will be in nats because we used the natural logarithm)
    entropy = -sum(logprobs .* exp(logprobs));
    
    % switch entropy from nats to bits
    entropy = entropy / log(2);
    
    
else    
    samples = input;
    
    % for a set of samples, we just delegate to the getEmpiricalModel which already does everything we need
    model = getEmpiricalModel(samples);
    entropy = model.entropy;
end


end