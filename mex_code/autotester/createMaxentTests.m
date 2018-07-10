% code that uses the old package for reference tests
function all_testdata = createMaxentTests()

load example15.mat
x = spikes15;

ncells = 15;
model_strings = {'indep','ising','kising','merp','ksync'};

input_models = {};
test_models = {};
all_testdata = {};
testnames = {};


% create models to test:

% regular models
for i = 1:numel(model_strings)
    model_input = maxent.createModel(ncells,model_strings{i});
    model_input = maxent.trainModel(model_input,x,'threshold',2);
    input_models{end+1} = model_input;
    test_models{end+1} = model_input; 
    testnames{end+1} = model_strings{i};
end

% test highorder model against ising model
correlations = cat(1,num2cell(nchoosek(1:ncells,1),2), ...
                num2cell(nchoosek(1:ncells,2),2));  % equivalent to the pairwise model
model_input = maxent.createModel(ncells,'ising');
model_input = maxent.trainModel(model_input,x,'threshold',2);
model_test = maxent.createModel(ncells,'highorder',correlations);
model_test = maxent.setFactors(model_test,maxent.getFactors(model_input));  % make them share the factors
model_test.z = model_input.z;
input_models{end+1} = model_input;
test_models{end+1} = model_test; 
testnames{end+1} = 'highorder';

% test composite model against kising
model_input = maxent.createModel(ncells,'kising');
model_input = maxent.trainModel(model_input,x,'threshold',2);
model_test1 = maxent.createModel(ncells,'ksync');
model_test2 = maxent.createModel(ncells,'ising');
model_test = maxent.createModel(ncells,'composite',{model_test1,model_test2});
model_test = maxent.setFactors(model_test,maxent.getFactors(model_input));  % make them share the factors
model_test.z = model_input.z;
input_models{end+1} = model_input;
test_models{end+1} = model_test; 
testnames{end+1} = 'composite';


% test sparse merp vs. regular merp
model_input = maxent.createModel(ncells,'merp');
model_input = maxent.trainModel(model_input,x,'threshold',2);
model_test = model_input;
model_input.type = 'merp';
model_test.type = 'merpsparse';
input_models{end+1} = model_input;
test_models{end+1} = model_test; 
testnames{end+1} = 'dense vs. sparse merp';


% tests for different models
for i = 1:numel(input_models)
    model_in = input_models{i};
    model_test = test_models{i};

    fprintf('*** Model_in: %s\n',model_in.type');
    
    % test marginals
    testdata.model = model_test;
    testdata.payload = ref_mexExhaustiveMarginals(model_in);
    testdata.command = 'data = maxent.getMarginals(testdata.model);';
    testdata.test = 'max(abs(testdata.payload-data))< epsilon';
    testdata.description = sprintf('marginals for %s',testnames{i});
    all_testdata{end+1} = testdata;
    
    % test gibbs
    params.randseed = 123;
    nsamples = 100;
    testdata.model = model_test;
    testdata.randseed = params.randseed;
    testdata.nsamples = nsamples;
    testdata.payload = ref_mexGibbsSampler(model_in,nsamples,0,params);
    testdata.command = 'xsynth = maxent.generateSamples(testdata.model,testdata.nsamples,''randseed'',testdata.randseed);';
    testdata.test = 'mean(xsynth(:) == testdata.payload(:))==1;';
    testdata.description = sprintf('sampler for %s',testnames{i});
    all_testdata{end+1} = testdata;

    % test empirical marginals
    testdata.model = model_test;
    testdata.x = x;
    testdata.payload = ref_mexEmpiricalMarginals(testdata.x,model_in);
    testdata.command = 'data = maxent.getEmpiricalMarginals(testdata.x,testdata.model);';
    testdata.test = 'max(abs(testdata.payload-data))< epsilon';
    testdata.description = sprintf('empirical marginals for %s',testnames{i});
    all_testdata{end+1} = testdata;

    % test getLogProbabilities
    testdata.model = model_test;
    testdata.x = x;
    testdata.payload = ref_mexLogProbability(testdata.x,model_in);
    testdata.command = 'data = maxent.getLogProbability(testdata.model,testdata.x);';
    testdata.test = 'max(abs(testdata.payload-data))< epsilon';
    testdata.description = sprintf('getlogprobability for %s',testnames{i});
    all_testdata{end+1} = testdata;
    
    % test getWeightedProbabilities
    testdata.model = model_test;
    testdata.x = x;
    testdata.p = exp(ref_mexLogProbability(testdata.x,model_in));
    testdata.payload = ref_mexEmpiricalMarginals(testdata.x,model_in,testdata.p);
    testdata.command = 'data = maxent.getWeightedMarginals(testdata.x,testdata.model,testdata.p);';
    testdata.test = 'max(abs(testdata.payload-data))< epsilon';
    testdata.description = sprintf('getWeightedMarginals for %s',testnames{i});
    all_testdata{end+1} = testdata;
    
    % test wang landau - except for indep where it is pointless and slow
    if (~strcmpi(model_test.type,'indep'))
        testdata.model = model_test;
        testdata.command = 'model_wang = rmfield(testdata.model,''z'');model_wang = maxent.wangLandau(model_wang);';
        testdata.test = 'abs(model_wang.z - testdata.model.z) < 5e-2';  % relatively lax threshold
        testdata.description = sprintf('wanglandau for %s',testnames{i});
        all_testdata{end+1} = testdata;
    end
    
    % test for trainModelExhaustive
    zero_factors = maxent.getFactors(model_test);
    zero_factors(:) = 0;
    model_untrained = maxent.setFactors(model_test,zero_factors);
    testdata.x = x;
    testdata.command = 'testdata.model = maxent.trainModelExhaustive(testdata.model,testdata.x,''silent'',true,''threshold'',2);';
    testdata.test = 'testMarginals(testdata.model,testdata.x,2.1);';
    testdata.description = sprintf('trainModelExhaustive for %s',testnames{i});
    all_testdata{end+1} = testdata;

    % test for trainModelMCMC
    zero_factors = maxent.getFactors(model_test);
    zero_factors(:) = 0;
    model_untrained = maxent.setFactors(model_test,zero_factors);
    testdata.x = x;
    testdata.command = 'testdata.model = maxent.trainModelMCMC(testdata.model,testdata.x,''silent'',true,''threshold'',2);';
    testdata.test = 'testMarginals(testdata.model,testdata.x,3);'; % less strict here because of MCMC
    testdata.description = sprintf('trainModelMCMC for %s',testnames{i});
    all_testdata{end+1} = testdata;


end


end