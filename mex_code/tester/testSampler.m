% tests the gibbs sampler on a specific model (compares to the reference)
function [out_results] =  testSampler(model,params,verbose)

EPSILON = 10e-7;

if nargin<3
    % print the results of the tests
    verbose = true;
end

bexact = true;
bapprox = true;


NSTEPS = 10000;


% test with 0 as starting point (inputed as a constant)
x_original = ref_gibbsSampler(0,NSTEPS,model,params);
x_new = mexGibbsSampler(model,NSTEPS,0,params);


fprintf('Gibbs / 0: ');
if (sum(x_original(:) ~= x_new(:)) == 0)
    printStatus('pass');
else
    bexact = false;
    bapprox = false;    
    printStatus('fail');
end

% test with a random starting point
x0 = rand(1,model.ncells);
x0 = uint32(x0 > 0.5);
x0_for_testing_orig = x0 * 1;  % force matlab to copy by value
x0_for_testing_new = x0 * 1;  % force matlab to copy by value
x_original = ref_gibbsSampler(x0_for_testing_orig,NSTEPS,model,params);
x0_original_result = x0_for_testing_orig;
x0_for_testing_new = x0;
x_new = mexGibbsSampler(model,NSTEPS,x0,params);
x0_new_result = x0_for_testing_new;
fprintf('Gibbs / random x0: ');
if (sum(x_original(:) ~= x_new(:)) == 0)
    printStatus('pass');
else
    bexact = false;
    bapprox = false;
    printStatus('fail');
end

fprintf('Gibbs / x0 modification: ');
if ((mean(x0_original_result==x0_new_result))==1)
    printStatus('pass');
else
    bexact = false;
    bapprox = false;
    printStatus('fail');
end


out_results.bexact = bexact;
out_results.bapprox = bapprox;


end