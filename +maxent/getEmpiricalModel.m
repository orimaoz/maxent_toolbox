% Returns an empirical distribution corresponding to a set of samples
%
% Usage:
%   model_out = getEmpiricalModel(samples)
%   model_out = getEmpiricalModel(samples,Name,Value,...)
%
%
% Arguments (mandatory):
%   samples - Set of samples, in the format (ncells x nsamples). 
%
% Optional arguments (in the form of Name,Value pairs):
%   min_count        - ignore all samples that appears less than this number of times (default 1)
%
% Output:
%   model_out - structure with several fields describing the empirical distrbution:
%       model_out.words  - unique codewords in the empirical distribution
%       model_out.counts - how many times each of these codewords repeated in the empirical dataset
%       model_out.logprobs - log probabilities of the codewords in the empirical distribution
%       model_out.logprobs - entropy of the empirical distribution, in bits
%

% Last update: Ori Maoz 09/11/2016
function model_out = getEmpiricalModel(samples,varargin)

DEFAULT_MIN_COUNT = 1;

p = inputParser;
addOptional(p,'min_count',DEFAULT_MIN_COUNT,@isnumeric);   % force exhaustive solver
p.parse(varargin{:});
min_count = p.Results.min_count;

[ncells,npatterns] = size(samples);

% find the different input words and their corresponding codewords
[unique_words I J] = unique(samples','rows','sorted');
counts = hist(J,numel(I));  % appearance count of each word

logprobs = log(counts / npatterns);

% keep only patterns for which we have reliable measurements
reliable_patterns = counts >= min_count;
logprobs = logprobs(reliable_patterns);
unique_words = unique_words(reliable_patterns,:);
counts = counts(reliable_patterns);

model_out.type = 'empirical';
model_out.ncells = ncells;
model_out.words = unique_words';
model_out.counts = counts;
model_out.logprobs = logprobs;
model_out.nsamples = npatterns;
model_out.entropy = -sum( (counts/npatterns) .* logprobs) / log(2);


end