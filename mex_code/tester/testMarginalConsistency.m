% Compares the mexExhaustiveMarginals function with mexGetMarginals to see if they are consistent with each other
function [out_results] =  testMarginalConsistency(model,verbose)

EPSILON = 10e-7;
NSTEPS = 10000;

if nargin<2
    % print the results of the tests
    verbose = true;
end

bexact = true;
bapprox = true;



% make a group of all patterns
npatterns = (2^model.ncells);
all_patterns = dec2bin(0:(npatterns-1))-'0';
all_patterns = uint32(all_patterns)';


% get some probabilities so we can do a weighted sum
all_logprobs = mexLogProbability(all_patterns,model);
all_logprobs = all_logprobs - logsumexp(all_logprobs(:));
all_probs = exp(all_logprobs);


% get the set of marginals and compare
m_all = mexEmpiricalMarginals(all_patterns,model,all_probs);
m_exhaustive = mexExhaustiveMarginals(model);
fprintf('Marginal consistency: ');
if (sum(m_all(:) ~= m_exhaustive(:)) == 0)
    printStatus('pass');        
else
    bexact = false;

    if (sum(abs(m_all-m_exhaustive)) < EPSILON)        
        printStatus('partial');        
    else
        bapprox = false;    
        printStatus('fail');        

    end

end


out_results.bexact = bexact;
out_results.bapprox = bapprox;


end