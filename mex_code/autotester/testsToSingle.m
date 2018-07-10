function tests = testsToSingle(tests)

for i = 1:numel(tests)
    if (strcmpi(tests{i}.model.type,'merp') || strcmpi(tests{i}.model.type,'merpsparse'))
        tests{i}.model.threshold = single(tests{i}.model.threshold);
        tests{i}.model.connections = single(tests{i}.model.connections);
        fprintf('changing test %d\n',i);
    end
    
end


end