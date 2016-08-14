

function timeAllModels



fprintf('Ising 20 cells\n');
load model_ising20
timeGibbs(model,100000);
timeMarginals(model,1000000);
timeLogProbability(model,1000000);
timeExhaustive(model);

fprintf('Ising 80 cells\n');
load model_ising80
timeGibbs(model,10000);
timeMarginals(model,100000);
timeLogProbability(model,100000);


fprintf('MERP 20 cells\n');
load model_merp20
timeGibbs(model,100000);
timeMarginals(model,100000);
timeLogProbability(model,1000000);
timeExhaustive(model);

fprintf('MERP 80 cells\n');
load model_merp80
timeGibbs(model,10000);
timeMarginals(model,10000);
timeLogProbability(model,10000);


fprintf('MERP sparse 20 cells\n');
load model_merpsparse20
timeGibbs(model,100000);
timeMarginals(model,100000);
timeLogProbability(model,1000000);
timeExhaustive(model);

fprintf('MERP sparse 80 cells\n');
load model_merpsparse80
timeGibbs(model,10000);
timeMarginals(model,10000);
timeLogProbability(model,10000);


fprintf('K-Ising 20 cells\n');
load model_kising20
timeGibbs(model,100000);
timeMarginals(model,1000000);
timeLogProbability(model,1000000);
timeExhaustive(model);

fprintf('K-Ising 80 cells\n');
load model_kising80
timeGibbs(model,10000);
timeMarginals(model,100000);
timeLogProbability(model,100000);


end