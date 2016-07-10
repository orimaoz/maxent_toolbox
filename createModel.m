% Initialize an empty maximum entropy model
%
% Usage:
%   model = createModel(ncells,model_string)
%   model = createModel(ncells,model_string,Name,Value,...)
%
% Arguments (mandatory):
%   ncells       - number of cells (dimensions) in the distribution
%   model_string - string denoting the type of models. Currently supported model types:
%                  indep  - independent model
%                  ising  - pairwise ising
%                  ksync  - k-synchrony
%                  kising - pairwise ising with k-synchrony
%                  merp   - MERP (maximum entropy based on random projections)
%
% Optional arguments (in the form of Name,Value pairs):
%   nprojections - number of random projections (for MERP model)
%   distribution - distribution the random projection values are drawn from (for MERP model)
%   distribution - distribution the random projection values are drawn from (for MERP model)
%   sparsity     - sparsity between 0 (completely sparse) and 1 (not sparse at all) - (for MERP model)
%   threshold    - relative threshold for random projection (for MERP model)
function model = createModel(ncells, model_string, varargin)

if nargin<2
    error('Usage: model = createModel(number_of_cells, model_type)')
end


% constants used for gradient descent in the models
STEP_SCALING_ISING = 1/10;
STEP_SCALING_KSYNC = 1;
STEP_SCALING_MERP = 1/50; 
STEP_SCALING_INDEP= 1; 


switch(lower(model_string))
    case {'indep'}   % independent model
        model.ncells = ncells;
        model.type = 'indep';
        model.factors = zeros(1,ncells);        
        model.step_scaling = ones(size(model.factors)) * STEP_SCALING_INDEP;
    case {'ising'}   % pairwise ME model
        model.ncells = ncells;
        model.type = 'ising';
        model.factors = zeros(1,ncells*(ncells+1)/2);
        model.step_scaling = ones(size(model.factors)) * STEP_SCALING_ISING;
    case {'kising'}  % pairwise ME model with synchrony constraint
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
    case {'merp'}    % MERP - maximum entropy with random projections
        
        % parse our optional arguments
        p = inputParser;
        addOptional(p,'nprojections',(ncells*(ncells+1)/2));  % default number of projections
        addOptional(p,'sparsity',5/ncells);                   % default sparsity 
        addOptional(p,'threshold',0.1);                       % default threshold 
        addOptional(p,'distribution','normal');               % default distribution of random projection values
        p.parse(varargin{:});
                
        model = createMERP(ncells,p.Results.nprojections,p.Results.distribution,p.Results.sparsity,p.Results.threshold);       
        model.step_scaling = ones(size(model.factors)) * STEP_SCALING_MERP;
    otherwise
        error('unsupported model string');
        
end


end