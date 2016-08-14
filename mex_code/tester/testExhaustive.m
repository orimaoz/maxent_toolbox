% tests the mexExhaustiveMarginals function on a specific model (compares to the reference)
function [out_results] =  testExhaustive(model,verbose)

EPSILON = 10e-7;
NSTEPS = 10000;

if nargin<2
    % print the results of the tests
    verbose = true;
end

bexact = true;
bapprox = true;



% get the set of marginals and compare
m_original = ref_mexExhaustiveMarginals(model);
m_new = mexExhaustiveMarginals(model);
fprintf('Exhaustive marginals: ');
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