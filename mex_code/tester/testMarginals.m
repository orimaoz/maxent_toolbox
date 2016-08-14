% tests the mexGetMarginals function on a specific model (compares to the reference)
function [out_results] =  testMarginals(model,verbose)

EPSILON = 10e-7;
NSTEPS = 10000;

if nargin<2
    % print the results of the tests
    verbose = true;
end

bexact = true;
bapprox = true;


% get some samples using one of the gibbs samplers
x = ref_gibbsSampler(0,NSTEPS,model);



% get the set of marginals and compare
m_original = ref_mexGetMarginals(x,model);
m_new = mexEmpiricalMarginals(x,model);
fprintf('Marginals / empirical: ');
if (sum(m_original(:) ~= m_new(:)) == 0)
    printStatus('pass');    
else
    bexact = false;

    if (sum(abs(m_original-m_new)) < EPSILON)        
        printStatus('partial');

    else
        bapprox = false;    
        printStatus('fail');
    end

end


% get some probabilities so we can do a weighted sum
lp_original = ref_mexLogProbability(x,model);
lp_original = lp_original - logsumexp(lp_original);
p_original = exp(lp_original);

m_original = ref_mexGetMarginals(x,model,p_original);
m_new = mexEmpiricalMarginals(x,model,p_original);
fprintf('Marginals / weighted: ');
if (sum(m_original(:) ~= m_new(:)) == 0)
    printStatus('pass');
else
    bexact = false;

    if (sum(abs(m_original-m_new)) < EPSILON)        
        printStatus('partial');
    else
        bapprox = false;    
        printStatus('fail');
    end

end


out_results.bexact = bexact;
out_results.bapprox = bapprox;


end