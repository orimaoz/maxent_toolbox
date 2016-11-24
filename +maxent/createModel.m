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
%                  composite- composite model combining the energy function of several other maximum entropy models.
%                             these should be provided in the 3rd argument as pre-created models in a cell array.
%
function model = createModel(ncells, model_string, varargin)

if nargin<2
    error('Usage: model = createModel(number_of_cells, model_type)')
end


% constants used for gradient descent in the models
STEP_SCALING_ISING = 1/10;
STEP_SCALING_KSYNC = 1;
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