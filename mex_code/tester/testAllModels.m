% test suite that compares the boltzmanntools mex function results to the "baseline" results to make sure that we 
% are not breaking anything while changing versions
function testAllModels


NCELLS = 15;
%RANDOM_SEED = 15;
RANDOM_SEED = randi([0 100000]);
NBURNIN = 0;

final_result_accurate = true;
final_result_approx = true;


params.burnin = NBURNIN;
params.randseed = RANDOM_SEED;



% test ising model
disp('*** TESTING ISING ***')
model = load('model_ising15');
collectResults(testSampler(model.model,params));
collectResults(testLogProbability(model.model));
collectResults(testMarginals(model.model));
collectResults(testExhaustive(model.model));
collectResults(testMarginalConsistency(model.model));
%collectResults(testWLS(model.model));


% test MERP
disp('*** TESTING MERP ***')
model = load('model_merp15');
model_merp = model;
collectResults(testSampler(model.model,params));
collectResults(testLogProbability(model.model));
collectResults(testMarginals(model.model));
collectResults(testExhaustive(model.model));
collectResults(testMarginalConsistency(model.model));
%collectResults(testWLS(model.model));

% test MERP
disp('*** TESTING SPARSE MERP ***')
model = load('model_merpsparse15');
collectResults(testSampler(model.model,params));
collectResults(testLogProbability(model.model));
collectResults(testMarginals(model.model));
collectResults(testExhaustive(model.model));
collectResults(testMarginalConsistency(model.model));
%collectResults(testWLS(model.model));

disp('*** COMPARING MERP WITH SPARSE MERP ***')
collectResults(testEquivalentModels(model.model,model_merp.model,params));

% test ksync
disp('*** TESTING K-synchrony ***')
model = load('model_ksync15');
collectResults(testSampler(model.model,params));
collectResults(testLogProbability(model.model));
collectResults(testMarginals(model.model));
collectResults(testExhaustive(model.model));
collectResults(testMarginalConsistency(model.model));
%collectResults(testWLS(model.model));

% test kising
disp('*** TESTING K-Ising ***')
model = load('model_kising15');
collectResults(testSampler(model.model,params));
collectResults(testLogProbability(model.model));
collectResults(testMarginals(model.model));
collectResults(testExhaustive(model.model));
collectResults(testMarginalConsistency(model.model));
%collectResults(testWLS(model.model));


disp('*** TESTING independent ***')
model = load('model_indep15');
collectResults(testSampler(model.model,params));
collectResults(testLogProbability(model.model));
collectResults(testMarginals(model.model));
collectResults(testExhaustive(model.model));
collectResults(testMarginalConsistency(model.model));



disp('*****************************');
fprintf('All tests: ');
if (final_result_accurate)
    printStatus('pass');
else
    if (final_result_approx)
        printStatus('partial');        
    else
        printStatus('fail');
    end        
end


    % helper function to collect the results (so that we don't miss a failure)
    function collectResults(results)
        final_result_accurate = final_result_accurate && results.bexact;         
        final_result_approx = final_result_approx && results.bapprox;         
    end


end