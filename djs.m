% Computes the Jensen-Shannon divergence between two distributions with a finite number of states.
% https://en.wikipedia.org/wiki/Jensen-Shannon_divergence
% The inputs are vectors of probabilities or log-probabilities (natural logarithm), the type of input is auto-detected.
% The probability distributions must be normalized (sum to 1), this is not checked by the function.
% 
% Usage: djs(p1,p2)
%
% Input:
%   p1 - Vector of probabilities for the first distribution (or log probabilities in natural logarithm).
%   p2 - Vector of probabilities for the second distribution (probabilities or natural log, same format as p1).
%
% Output:
%   Jensen-Shannon divergence between the distributions: DJS(p1||p2). 
function out = djs(p1,p2)

if (min(p1) < 0)
    % we were given log probabilities, we must convert them to regular probabilities before we can sum them
    p1 = exp(p1);
    p2 = exp(p2);
end

% DJS is just DKL of the two distributions to a middle distribution
m = (p1+p2)/2;  
out = (dkl(p1,m) + dkl(p2,m))/2;

end