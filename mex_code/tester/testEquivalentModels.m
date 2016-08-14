% Tests two models that should be the same (i.e. sparse/nonsparse versions) to see that they are consistent
function [out_results] =  testEquivalentModels(model1,model2,params)

EPSILON = 10e-7;


bexact = true;
bapprox = true;
NSTEPS = 1000;


% get the set of marginals and compare
m1 = mexExhaustiveMarginals(model1);
m2 = mexExhaustiveMarginals(model2);
fprintf('Exhaustive marginals: ');
if (sum(m1(:) ~= m2(:)) == 0)
    printStatus('pass');    
else
    bexact = false;

    if (sum(abs(m1-m2)) < EPSILON)        
        printStatus('partial');
    else
        bapprox = false;    
        printStatus('fail');
    end

end


% test gibbs sampler with 0 as starting point (inputed as a constant)
x1 = mexGibbsSampler(model1,NSTEPS,0,params);
x2 = mexGibbsSampler(model2,NSTEPS,0,params);


fprintf('Gibbs / 0: ');
if (sum(x1(:) ~= x2(:)) == 0)
    printStatus('pass');
else
    bexact = false;
    bapprox = false;    
    printStatus('fail');
end


out_results.bexact = bexact;
out_results.bapprox = bapprox;


end