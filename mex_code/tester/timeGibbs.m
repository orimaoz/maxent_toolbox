function [old_time, new_time] = timeGibbs(model,nsteps)

% make the starting point random but common for both
params.randseed = randi([0 100000]);
params.burnin = 0;

a = tic;
x = ref_gibbsSampler(0,nsteps,model,params);
old_time = toc(a);


a = tic;
x = gibbsSampler(0,nsteps,model,params);
new_time = toc(a);


fprintf('Gibbs     ');
printTiming(old_time,new_time);

end