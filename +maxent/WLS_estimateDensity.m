% Implementation of Wang Landau sampling.
% invoked by wangLandau() and typically should not be invoked directly.
%
% 
% based on:
% "Efficient, Multiple-Range Random Walk Algorithm to Calculate the Density of States",
% Wang and Landau (2001), Physical Review Letters.
%
% input:
% model - model to evaluate
% nstates - number of states in histogram
% niterations - how many iterations to perform (final value of f is 2^(-iterations))
% separation - separation between sampled values in MCMC walk
%
% output:
% e - histogram centers for energies
% g - densities of states
%
% Last update: Ori Maoz, July 2016
function [e g] = WLS_estimateDensity(model,varargin)

p = inputParser;
addRequired(p,'binsize'); % energy bin size
addRequired(p,'depth'); % depth of MCMC interations
addRequired(p,'separation'); % separation between samples
addRequired(p,'savefile');                 % save file for re-entrant code
addRequired(p,'save_delay');    % time between re-entry file saves
addOptional(p,'silent',false,@islogical);    % silent mode - don't print anything


p.parse(varargin{:});
delta = p.Results.binsize;
niterations = p.Results.depth;
separation = p.Results.separation;
reentry_filename = p.Results.savefile;
time_between_saves = p.Results.save_delay;
silent = p.Results.silent;


FLATNESS_CRITEREA = 0.9;


ncells = model.ncells;

factors = maxent.getFactors(model);





max_energy = sum(factors(factors>0));
min_energy = sum(factors(factors<0));
e = (min_energy:delta:(max_energy+delta));





% make sure that the reentry filename is well-structured
last_save = tic;
if ~isempty(reentry_filename)
    [pathstr,name,~] = fileparts(reentry_filename);
    reentry_filename  = fullfile(pathstr,[name '.mat']);
    
    internal_print('Using save file %s for re-entry\n',reentry_filename);
end



%if the re-entry file exists, load running status from it
if exist(reentry_filename, 'file') == 2
    load(reentry_filename);
    
    % reset the time of the last save so we won't re-save immediately
    last_save = tic;
    
else
    
    
    % compute the size of each bin and spread the energies around with fixed bins.
    persistent_state.energies.bins = double(min_energy:delta:max_energy);

    
    % apriori density of states set to 1
    persistent_state.energies.densities = zeros(size(persistent_state.energies.bins));
    persistent_state.energies.histogram = zeros(size(persistent_state.energies.bins));

    persistent_state.energies.update_factor = 1;    
    
    
    % Starting point for the random walk
    persistent_state.x = uint32(zeros(1,ncells));
    
    persistent_state.layer = 1;
end




fprintf(['[' repmat('.',[1 (niterations)]) ']']);
% for layer = 1:niterations
while persistent_state.layer <= niterations
    fprintf(repmat('\b',[1 (niterations+1)]));
    fprintf([repmat('.',[1 (persistent_state.layer)]) repmat(' ',[1 (niterations-persistent_state.layer)]) ']']);
         
    % reset the histogram
    persistent_state.energies.histogram = zeros(size(persistent_state.energies.bins));
    persistent_state.flat_histogram = false;

    while (~persistent_state.flat_histogram)

        
        
        mexWangLandauSampler(persistent_state.x,500000,persistent_state.energies,model,separation);
        
        % normalize the histgram by the smallest nonzero entry to retain numeric precision
        nonzero_entries = (persistent_state.energies.densities>0);

        
        % check if the histogram is flat
        minH = min(persistent_state.energies.histogram(nonzero_entries));
        meanH = mean(persistent_state.energies.histogram(nonzero_entries));
        
        %disp(sprintf('Histogram minimum/mean: %.03f/%.03f',minH,meanH));
        if (minH > meanH * FLATNESS_CRITEREA)
 %          disp(sprintf('Finished with f=%f  (logf = %g)',exp(energies.update_factor),energies.update_factor));
            persistent_state.flat_histogram = true; 
        end

        % check if we should save re-entry data
        if ~isempty(reentry_filename)
            time_from_last_save = toc(last_save);
            if (time_from_last_save > time_between_saves)
                % save the variables necessary to resume operation if we are killed and re-started with the same input.
                save(reentry_filename,'persistent_state');

                last_save = tic;
            end
        end

        
    end

    % next iteration - smaller modification factor
    persistent_state.energies.update_factor = persistent_state.energies.update_factor/2;
    
    
    persistent_state.layer = persistent_state.layer + 1;
end % f_layers
fprintf('\n');

    

    % remove all the parts of the histogram that were never visited
    visited_states = (persistent_state.energies.densities>0);
    g = persistent_state.energies.densities(visited_states);
    e = persistent_state.energies.bins(visited_states);
    
    function internal_print(varargin)
        if (~silent)
            disp(sprintf(varargin{:}));
        end        
    end

end
