% Initialize an empty maximum entropy model
%
% Usage:
%   model = createModel(ncells,'model_string','param1',value1,'param2',value2);
%
% Supported model strings:
%   indep  - independent model
%   ising  - pairwise ising
%   ksync  - k-synchrony
%   kising - pairwise ising with k-synchrony
%   merp   - MERP (maximum entropy based on random projections)
%
% The optional value pairs can be either given as arguments or mashed together in a structure.
% they contain model-specific information
%   
%   merp:
%   nprojections - number of random projections
%   distribution - distribution the random projection values are drawn from
%   sparsity - sparsity between 0 (completely sparse) and 1 (not sparse at all)
%   threshold - threshold for firing
function model = createModel(ncells, model_string, varargin)

% constants used for gradient descent in the models
STEP_SCALING_ISING = 1/10;
STEP_SCALING_KSYNC = 1;
STEP_SCALING_MERP = 1/50; 
STEP_SCALING_INDEP= 1; 



switch(lower(model_string))
    case {'ising'}
        model.ncells = ncells;
        model.type = 'ising';
        model.factors = zeros(1,ncells*(ncells+1)/2);
        model.step_scaling = ones(size(model.factors)) * STEP_SCALING_ISING;
    case {'kising'}
        model.ncells = ncells;
        model.type = 'kising';
        model.factors = zeros(1,(ncells+2)*(ncells+1)/2);
        model.step_scaling = ones(size(model.factors)) * STEP_SCALING_ISING;
        model.step_scaling(1:(ncells+1)) = STEP_SCALING_KSYNC;
    case {'ksync'}
        model.ncells = ncells;
        model.type = 'ksync';
        model.factors = zeros(1,(ncells+1));
        model.step_scaling = ones(size(model.factors)) * STEP_SCALING_KSYNC;
    case {'merp'}
        
        % parse our optional arguments
        p = inputParser;
        addOptional(p,'nprojections',(ncells*(ncells+1)/2));  % default number of projections
        addOptional(p,'sparsity',3/ncells);                   % default sparsity 
        addOptional(p,'threshold',0.1);                       % default threshold 
        addOptional(p,'distribution','normal');               % default distribution of random projection values
        p.parse(varargin{:});
                
        model = createMERP(ncells,p.Results.nprojections,p.Results.distribution,p.Results.sparsity,p.Results.threshold);       
        model.step_scaling = ones(size(model.factors)) * STEP_SCALING_MERP;
    case {'indep'}
        model.ncells = ncells;
        model.type = 'indep';
        model.factors = zeros(1,ncells);        
        model.step_scaling = ones(size(model.factors)) * STEP_SCALING_INDEP;
    otherwise
        error('unsupported model string');
        
end


end