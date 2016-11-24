% Internal function for creating a MERP model.
% this function should not be invoked directly, instead call the "createModel" function.
% parameters:
%   ncells - number of input cells
%   nfactors - number of projections
%   distributions - distribution to draw the weights from
%   sparsity - sparsity of connections, from 1 (all connections exist) to 0 (no connections at all)
%   threshold - firing threshold, per connected cell.
%   exact - if this is set to true, each cell will have exactly (ncells*sparsity) connections, otherwise each
%           connection is pruned with a probability of (1-sparsity)
%
function net = createMERP(ncells,nfactors,distribution,sparsity,threshold,exact)

if nargin<6
    exact = false;
end

if (sparsity <=0 || sparsity > 1)
    error('MERP sparsity must be in the range (0,1]');
end

% we can draw the factors from a variety of distributions
switch lower(distribution)
    case 'powerlaw'
        pd = makedist('Gamma',3,2);
        bias = -6;
        net.connections = (pd.random(nfactors,ncells)) + bias;      
    case 'lognormal'
        pd = makedist('LogNormal');
        net.connections = (pd.random(nfactors,ncells));      
    case 'uniform'
        % uniform non-zero distribution
        pd = makedist('uniform');
        net.connections = (pd.random(nfactors,ncells)) - 0.2;      
    case 'binary'
        % all active synapses are the same strength
        pd = makedist('multinomial','Probabilities',[0.2, 0.8]);
        net.connections = (pd.random(nfactors,ncells) - 1.5)*2;
    case 'normal'
        % normal distribution which can also have negative synapses
        pd = makedist('normal',1,1);
        net.connections = (pd.random(nfactors,ncells));      
    otherwise
        error('Unknown connection distribution inputed');
end
       
% randomly draw synapses according to the supplied distribution
%net.connections = (pd.random(nfactors,ncells)) + bias;        
net.factors = zeros(1,nfactors);
net.ncells = ncells;    

% threshold is relative to how many cells participate and how strong the connetions are
net.threshold = ones(1,nfactors) * threshold * ncells * sparsity * mean(net.connections(:));    
net.z = 0;

% if the model is sparse then we take advantage of the sparsity. If it is not sparse enough, it is better
% to work in a non-sparse manner. We choose this according to a ballpark estimate (TODO: measure speeds later)
if sparsity < 0.4   
    net.type = 'merpsparse';
else
    net.type = 'merp';
end


if (exact)
    % this version leaves a fixed number of synapses from each cell:
    %nremaining = ceil((1-sparsity)*ncells);
    nremove = floor((1-sparsity)*ncells);
    for i = 1:nfactors
        p = randperm(ncells);
        net.connections(i,p(1:nremove)) = 0;       
    end
else
    % apply sparsity - remove some of the connections
    % this version only gives a removal probability for each syanpse:
    if (sparsity<1)
        remove_mask = rand(size(net.connections));
        net.connections(remove_mask > sparsity) = 0;
    end % sparsity
end 

        
end