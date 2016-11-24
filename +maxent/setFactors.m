% sets the factors of a model from an array.
%
% Usage:
%   model = createModel(model)
%
% returns:
%   model structure with updated factors
function model = setFactors(model,factors)
    
    if (strcmpi(model.type,'composite'))
        % if it's a composite model, recursively iterate over all the submodels
        factors_idx = 1;
        for i = 1:numel(model.innermodels)
            innermodel = model.innermodels{i};            
            nfactors_inner = maxent.getNumFactors(innermodel);
            
            model.innermodels{i} = maxent.setFactors(model.innermodels{i},factors(factors_idx:(factors_idx+nfactors_inner-1)));
            
            factors_idx = factors_idx + nfactors_inner;
            
        end
        
    else
        % if it's a simple model, just set the factors
        model.factors = factors;
    end

end