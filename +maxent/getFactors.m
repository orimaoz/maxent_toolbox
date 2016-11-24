% returns the factors of a model as an array
%
% Usage:
%   model = getFactors(model)
function factors = getFactors(model)
    
    if (strcmpi(model.type,'composite'))
        % if it's a composite model, recursively iterate over all the submodels
        factors = [];
        for i = 1:numel(model.innermodels)
            factors = [factors, maxent.getFactors(model.innermodels{i})];
        end
        
    else
        % if it's a simple model, just return the factors
        factors = model.factors;
    end

end