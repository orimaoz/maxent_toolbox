% prints timing results for functions
function printTiming(old_time, new_time)

threshold_fast = 1.2;
threshold_slow = 1.1;

fprintf('Original: %f ',old_time);
fprintf('New : %f    ',new_time);

if ((old_time / new_time) > threshold_fast)
    cprintf([0,0.6,0],'(FAST)');    
end

if ((new_time / old_time) > threshold_slow)   
    cprintf([0.6,0,0],'(SLOW)');    
end

fprintf('\n');

end