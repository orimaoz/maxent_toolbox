% tests the wang landau sampler on a specific model (compares to the reference)
function [out_results] =  testWLS(model,verbose)

EPSILON = 10e-7;

if nargin<3
    % print the results of the tests
    verbose = true;
end

bexact = true;
bapprox = true;


NSTEPS = 10000;
ENERGY_DELTA = 0.1;


% set up parameters for wang landau sampling
max_energy = -sum(model.factors(model.factors<0));
min_energy = -sum(model.factors(model.factors>0));
e = (min_energy:ENERGY_DELTA:(max_energy+ENERGY_DELTA));
energies.bins = double(min_energy:ENERGY_DELTA:max_energy);
energies.densities = zeros(size(energies.bins));
energies.histogram = zeros(size(energies.bins));
energies.update_factor = 1;

% Starting point for the random walk
x = uint32(zeros(1,model.ncells));

% test with a random starting point
x0 = rand(1,model.ncells);
x0 = uint32(x0 > 0.5);
x0_for_testing_orig = x0 * 1;  % force matlab to copy by value
x0_for_testing_new = x0 * 1;  % force matlab to copy by value
ref_energies = energies;
ref_energies.densities = ref_energies.densities * 1;
ref_energies.histogram = ref_energies.histogram * 1;
ref_wangLandauSampler(x0_for_testing_orig,NSTEPS,ref_energies,model);


x0_original_result = x0_for_testing_orig;
new_energies = energies;
new_energies.densities = new_energies.densities * 1;
new_energies.histogram = new_energies.histogram * 1;
wangLandauSampler(x0_for_testing_new,NSTEPS,new_energies,model);
x0_new_result = x0_for_testing_new;
fprintf('Wanglandau / densities: ');
if (sum(new_energies.densities ~= ref_energies.densities) == 0)
    %disp('Wanglandau / densities: pass');
    printStatus('pass');            
else
    bexact = false;
    bapprox = false;
    %disp('Wanglandau / densities: FAIL');    
    printStatus('fail');        

end
fprintf('Wanglandau / histogram: ');
if (sum(new_energies.histogram ~= ref_energies.histogram) == 0)
    %disp('Wanglandau / histogram: pass');
    printStatus('pass');            
else
    bexact = false;
    bapprox = false;
    %disp('Wanglandau / histogram: FAIL');    
    printStatus('fail');            
end

fprintf('Wanglandau / x0 modification: ');
if ((mean(x0_original_result==x0_new_result))==1)
    %disp('Wanglandau / x0 modification: pass');
    printStatus('pass');            
else
    bexact = false;
    bapprox = false;
    %disp('Wanglandau / x0 modification: FAIL');        
    printStatus('fail');            
end


out_results.bexact = bexact;
out_results.bapprox = bapprox;


end