function [old_time, new_time] = timeMarginals(model,nsteps)


% create some data
x = rand(model.ncells, nsteps);
x = uint32(x > 0.5);


a = tic;
m = ref_mexGetMarginals(x,model);
old_time = toc(a);


a = tic;
m = mexGetMarginals(x,model);
new_time = toc(a);


fprintf('Marginals ');
printTiming(old_time,new_time);
end