% Normalizes a probabilistic model.
% If the model is of less than 30 dimensions, it is normalized exactly by summing over all its states.
% Otherwise, it is normalized approximately by calling the wangLandau function.
%
% Usage:
%   entropy = normalizeModel(model)
%
% Arguments (mandatory):
%   model       - Maximum entropy model as returned by trainModel()
%
% Output:
%   model_out   - same model after setting its 'z' (normalization) value so that it is properly normalized
%
% Last update: Ori Maoz 08/01/2018
function model_out = normalizeModel(model,varargin)


    if model.ncells < 30       
        % exhaustive (exact) computation
        
        % we will just shortcut and use the built-in mex implementation that also computes the marginals
        [marginals,z] = mexExhaustiveMarginals(model);
        model_out = model;
        model_out.z = -log(z);    
                       
    else
        % approximation based on the Wang-Landau implementation, using the default arguments or whatever arguments were
        % passed into this function
        model_out = maxent.wangLandau(model,varargin{:});
    end
    


end
