% Computes the Kullback-Leibler divergence between two distributions with a finite number of states.
% https://en.wikipedia.org/wiki/Kullback-Leibler_divergence
% The inputs are vectors of probabilities or log-probabilities (natural logarithm), the type of input is auto-detected.
% The probability distributions must be normalized (sum to 1), this is not checked by the function.
% 
% Usage: dkl(p1,p2)
%
% Input:
%   p1 - Vector of probabilities for the first distribution (or log probabilities in natural logarithm).
%   p2 - Vector of probabilities for the second distribution (probabilities or natural log, same format as p1).
%
% Output:
%   Kullback-Leibler divergence between the distributions: DKL(p1||p2). 

function out = dkl(p1,p2)

if (min(p1) < 0)
    % we were given log probabilities.
    
    % we assume that the distributions were given in the natural logarithm, so first convert them to the
    % base-2 logarithm so that the result will be in bits
    lp1 = p1(:) / log(2);
    lp2 = p2(:) / log(2);
    
    % now compute DKL
    ratios = 2.^(lp1) .* (lp1-lp2);
    out = nansum(ratios);
else
    % we were given probabilities
    p1 = p1(:);
    p2 = p2(:);
    
    ratios = p1.*log2(p1) - p1 .* log2(p2);
    out = nansum(ratios);
end


end