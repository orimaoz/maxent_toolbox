function all_pass = runMaxentTests(all_testdata)

% running tests

epsilon = 1e-10;
all_pass = true;
for i = 1:numel(all_testdata)
    testdata = all_testdata{i};
    fprintf('%s: ',testdata.description);
    eval(testdata.command);
    bresult = eval(testdata.test);
    if bresult
        fprintf('pass\n');
    else
        fprintf('FAIL\n');        
    end
    
    all_pass = all_pass && bresult;
end

if all_pass
    fprintf('all tests: pass\n');
else
    fprintf('all tests: FAIL\n');        
end

end