% Initialize an empty maximum entropy model
%
% Usage:
%   model = createModel(ncells,model_string)
%   model = createModel(ncells,'highorder',correlations)
%   model = createModel(ncells,'composite',inner_models)
%   model = createModel(ncells,model_string,Name,Value,...)
%
% Arguments (mandatory):
%   ncells       - number of cells (dimensions) in the distribution
%   model_string - string denoting the type of models. Currently supported model types:
%                  indep    - independent model
%                  pairwise - pairwise maximum entropy
%                  ksync    - k-synchrony
%                  kising   - pairwise ME with k-synchrony
%                  highorder- set of custom user-supplied interaction, supplied as third argument in the form of a cell
%                             array of arrays, for example: {[1 3],[3 5 7],[2 3]}
%                  rp       - RP (Random Projection) model
%                  composite- composite model combining the energy function of several other maximum entropy models.
%                             these should be provided in the 3rd argument as pre-created models in a cell array.
%
% Optional arguments (in the form of Name,Value pairs):
%   nprojections - number of random projections (for RP model)
%   distribution - distribution the random projection values are drawn from (for RP model)
%   indegree     - average number of nonzero projection elements for RP model
%   threshold    - relative threshold for random projection (for RP model)
%   interactions - cell array of interacting groups for high-order maxent model (for "custom" model)
function model = createModel(ncells, model_string, varargin)

if nargin<2
    error('Usage: model = createModel(number_of_cells, model_type)')
end


% constants used for gradient descent in the models
STEP_SCALING_ISING = 1/10;
STEP_SCALING_KSYNC = 1;
STEP_SCALING_RP = 1/50; 
STEP_SCALING_INDEP= 1; 


switch(lower(model_string))
    case {'indep'}   % independent model
        model.ncells = ncells;
        model.type = 'indep';
        model.factors = zeros(1,ncells);        
        model.step_scaling = ones(size(model.factors)) * STEP_SCALING_INDEP;
    case {'ising','pairwise'}   % pairwise ME model
        model.ncells = ncells;
        model.type = 'ising';
        model.factors = zeros(1,ncells*(ncells+1)/2);
        model.step_scaling = ones(size(model.factors)) * STEP_SCALING_ISING;
    case {'kising','kpairwise'}  % pairwise ME model with synchrony constraint
        model.ncells = ncells;
        model.type = 'kising';
        model.factors = zeros(1,(ncells+2)*(ncells+1)/2);
        model.step_scaling = ones(size(model.factors)) * STEP_SCALING_ISING;
        model.step_scaling(1:(ncells+1)) = STEP_SCALING_KSYNC;
    case {'ksync'}   % ME model with synchrony constraints
        model.ncells = ncells;
        model.type = 'ksync';
        model.factors = zeros(1,(ncells+1));
        model.step_scaling = ones(size(model.factors)) * STEP_SCALING_KSYNC;
    case {'highorder'} % User-supplied list of high-order interactions
        % get required argument - user-specified interactions
        p = inputParser;
        addRequired(p,'interactions',@iscell);  % default number of projections
        p.parse(varargin{:});
        interactions = p.Results.interactions;
        num_constraints = numel(interactions);
                
        % use a trick - implement custom interactions as a specific projection matrix of RP
        model.ncells = ncells;
        model.type = 'highorder';
        model.factors = zeros(1,num_constraints);
        model.step_scaling = ones(size(model.factors)) * STEP_SCALING_ISING;
        model.factorMatrix = single(zeros(num_constraints,ncells));
        model.threshold = single(zeros(1,num_constraints));
        
        for i = 1:num_constraints
            model.factorMatrix(i,interactions{i}) = 1;
            model.threshold(i) = numel(interactions{i}) - 0.5;
        end
        
    case {'rp','RP','merp'}    % RP - Random Projection model
        
        % parse our optional arguments
        p = inputParser;
        addOptional(p,'nprojections',(ncells*(ncells+1)/2));  % default number of projections
        addOptional(p,'sparsity',nan);                   % default sparsity 
        addOptional(p,'threshold',0.1);                       % default threshold 
        addOptional(p,'indegree',5);
        p.parse(varargin{:});
            
        nprojections = p.Results.nprojections;
        threshold = p.Results.threshold;
        indegree = p.Results.indegree;
        
        % sparsity of the model (average number of nonzero elements per input dimension) can be specified either
        % explicitly or implicitly (by choosing an average indegree)
        sparsity = indegree / ncells;
        if ~isnan(p.Results.sparsity)
            sparsity = p.Results.sparsity;
        end

        % draw synaptic values from a normal distribution
        pd = makedist('normal',single(1),single(1));
        model.A = (pd.random(nprojections,ncells));              
                
        % keep only some of the synapses to be nonzero according to the chosen sparsity
        if (sparsity<1)
            remove_mask = rand(size(model.A));
            model.A(remove_mask > sparsity) = 0;
        end % sparsity

        model.sparse = (sparsity < 0.4);
        
        model.type = 'rp';
        model.factors = zeros(1,nprojections);
        model.ncells = ncells;    
        
        % threshold is relative to how many cells participate 
        model.threshold = single(ones(1,nprojections) * threshold * ncells * sparsity);    
        %model.threshold = single(ones(1,nfactors) * threshold * ncells * sparsity * mean(net.connections(:)));                   
        model.z = 0;
                       
        model.step_scaling = ones(size(model.factors)) * STEP_SCALING_RP;
    
    case {'composite'}
        % it's a composite model made up of several smaller models
        p = inputParser;
        addRequired(p,'innermodels',@iscell);  % contained models
        p.parse(varargin{:});
        innermodels = p.Results.innermodels;
        
        if isempty(innermodels)
            error('list of contained models cannot be empty');
        end
        
        model.ncells = ncells;        
        model.type = 'composite';
        
        % check the inner models one by one to see that they are of the correct dimension
        for i = 1:numel(innermodels)
            if (innermodels{i}.ncells ~= ncells)
                error('Inner model has different number of dimensions than composite model');
            end
            
            model.innermodels{i} = innermodels{i};
            
        end
        
    
    otherwise
        error('unsupported model string');
        
end


end