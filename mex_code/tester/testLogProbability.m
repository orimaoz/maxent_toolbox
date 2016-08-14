% tests the mexLogProbability function on a specific model (compares to the reference)
function [out_results] =  testLogProbability(model,verbose)

EPSILON = 10e-7;
NSTEPS = 10000;
DEFAULT_Z_VALUE = 2;

if nargin<2
    % print the results of the tests
    verbose = true;
end

bexact = true;
bapprox = true;


% get some samples using one of the gibbs samplers
x = ref_gibbsSampler(0,NSTEPS,model);


% if the model is not normalized, force some normalization just so we can make sure it is working as expected
if isfield(model,'z')
    if model.z ~= 0
       model.z = DEFAULT_Z_VALUE; 
    end
else
    model.z = DEFAULT_Z_VALUE;
end


% test with 0 as starting point (inputed as a constant)
p_original = ref_mexLogProbability(x,model);
p_new = mexLogProbability(x,model);
fprintf('LogProbability / normalized: ');
if (sum(p_original(:) ~= p_new(:)) == 0)
    printStatus('pass');
    
else
    bexact = false;

    if (sum(abs(p_original-p_new)) < EPSILON)        
        printStatus('partial');

    else
        bapprox = false;    
        printStatus('fail');
    end

end





% test with 0 as starting point (inputed as a constant)
model = rmfield(model,'z');
p_original = ref_mexLogProbability(x,model);
p_new = mexLogProbability(x,model);
fprintf('LogProbability / unnormalized: ');
if (sum(p_original(:) ~= p_new(:)) == 0)
    printStatus('pass');

else
    bexact = false;

    if (sum(abs(p_original-p_new)) < EPSILON)        
        printStatus('partial');
    else
        bapprox = false;    
        printStatus('fail');
    end

end


out_results.bexact = bexact;
out_results.bapprox = bapprox;


end