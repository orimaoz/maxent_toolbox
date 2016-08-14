function [old_time, new_time] = timeLogProbability(model,nsteps)


% create some data
x = rand(model.ncells, nsteps);
x = uint32(x > 0.5);


a = tic;
m = ref_mexLogProbability(x,model);
old_time = toc(a);


a = tic;
m = mexLogProbability(x,model);
new_time = toc(a);


fprintf('Log prob  ');
printTiming(old_time,new_time);
end