% Returns an explicit representation of a probability distribution as a vector of probabilities.
% This function will return an error for large (n >= 30) distributions since they typically cannot fit in memory.
%
% Usage:
%   probabilities = getExplicitDistribution(model)
%
% Arguments (mandatory):
%   model       - Maximum entropy model as returned by trainModel()
%
% Output:
%   probabilities - vector of probabilities for all system states starting from 00000, 00001.... up to 11111
%
% Last update: Ori Maoz 28/01/2018
function probabilities = getExplicitDistribution(model)
    
    if model.ncells < 30       
        
        % get a vector of all possible states
        npatterns = (2^model.ncells); 
        all_states = logical(de2bi(0:(npatterns-1)))';

        % get the probability for each
        all_log_probabilities = maxent.getLogProbability(model,all_states);
                
        probabilities = exp(all_log_probabilities);
        
        % check if we need to normalize the distribution
        prob_sum = sum(probabilities);
        if (abs(prob_sum-1)>eps)
            probabilities = probabilities / prob_sum;
        end
        
                       
    else
        error('Cannot work on n >= 30, will take up too much memory!');
    end
    


end
