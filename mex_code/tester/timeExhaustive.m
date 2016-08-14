function [old_time, new_time] = timeExhaustive(model)


a = tic;
m = ref_mexExhaustiveMarginals(model);
old_time = toc(a);


a = tic;
m = mexExhaustiveMarginals(model);
new_time = toc(a);


fprintf('Exhastive ');
printTiming(old_time,new_time);


end