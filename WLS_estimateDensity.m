% Implementation of Wang Landau sampling.
% invoked by wangLandau() and typically should not be invoked directly.
%
% 
% based on:
% "Efficient, Multiple-Range Random Walk Algorithm to Calculate the Density of States",
% Wang and Landau (2001), Physical Review Letters.
%
% input:
% model - ising model to evaluate
% nstates - number of states in histogram
% niterations - how many iterations to perform (final value of f is 2^(-iterations))
% separation - separation between sampled values in MCMC walk
%
% output:
% e - histogram centers for energies
% g - densities of states
%
% Last update: Ori Maoz, July 2016
function [e g] = WLS_estimateDensity(model,delta,niterations,separation)


FLATNESS_CRITEREA = 0.9;

if nargin<3
    niterations = 20;
end


ncells = model.ncells;


max_energy = sum(model.factors(model.factors>0));
min_energy = sum(model.factors(model.factors<0));
e = (min_energy:delta:(max_energy+delta));


% compute the size of each bin and spread the energies around with fixed bins.
energies.bins = double(min_energy:delta:max_energy);


% apriori density of states set to 1
energies.densities = zeros(size(energies.bins));
energies.histogram = zeros(size(energies.bins));

% initial modification factor
energies.update_factor = 1;


% Starting point for the random walk
x = uint32(zeros(1,ncells));

fprintf(['[' repmat('.',[1 (niterations)]) ']']);
for layer = 1:niterations
    fprintf(repmat('\b',[1 (niterations+1)]));
    fprintf([repmat('.',[1 (layer)]) repmat(' ',[1 (niterations-layer)]) ']']);
         
    % reset the histogram
    energies.histogram = zeros(size(energies.bins));
    flat_histogram = false;

    while (~flat_histogram)

        
        
        mexWangLandauSampler(x,500000,energies,model,separation);
        
        % normalize the histgram by the smallest nonzero entry to retain numeric precision
        nonzero_entries = (energies.densities>0);

        
        % check if the histogram is flat
        minH = min(energies.histogram(nonzero_entries));
        meanH = mean(energies.histogram(nonzero_entries));
        
        %disp(sprintf('Histogram minimum/mean: %.03f/%.03f',minH,meanH));
        if (minH > meanH * FLATNESS_CRITEREA)
 %          disp(sprintf('Finished with f=%f  (logf = %g)',exp(energies.update_factor),energies.update_factor));
           flat_histogram = true; 
        end

    end

    % next iteration - smaller modification factor
    energies.update_factor = energies.update_factor/2;
end % f_layers
fprintf('\n');

    

    % remove all the parts of the histogram that were never visited
    visited_states = (energies.densities>0);
    g = energies.densities(visited_states);
    e = energies.bins(visited_states);

end
