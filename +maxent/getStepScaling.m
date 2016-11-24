% returns the step scaling for a model as an array
%
% Usage:
%   model = getFactors(model)
%
% returns: step scaling if it exists, or nan if not
function scaling = getStepScaling(model)
    
    if (strcmpi(model.type,'composite'))
        % if it's a composite model, recursively iterate over all the submodels
        scaling = [];
        for i = 1:numel(model.innermodels)
            scaling = [scaling, maxent.getStepScaling(model.innermodels{i})];
        end
        
    else
        % if it's a simple model, just return the factors
        if isfield(model,'step_scaling')
            scaling = model.step_scaling;
        else
            scaling = nan;
        end
    end

end