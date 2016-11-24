% returns the (total) number of factors in a model
%
% Usage:
%   model = getNumFactors(model)
function nfactors = getNumFactors(model)
    
    if (strcmpi(model.type,'composite'))
        % if it's a composite model, recursively iterate over all the submodels
        nfactors = 0;
        for i = 1:numel(model.innermodels)
            nfactors = nfactors + maxent.getNumFactors(model.innermodels{i});
        end
        
    else
        % if it's a simple model, just return the factors
        nfactors = numel(model.factors);
    end

end