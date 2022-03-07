function [Os, X] = binary_constraints_multipartite(n, K)
% Inputs: Number of binary variables n in the system, and order K of
% marginal transition probabilities to constrain.
%
% Outputs: 2^n x n x m  array of constrained operators Os, where
% Os(:,:,i) is the i^th constraint; and 2^n x n matrix of states X, where
% X(i,:) is the i^th state. We note that Os(i,j,:) represents the
% constraints on variable j flipping in state i.

% Matrix of states:
X = double(dec2bin(0:(2^n - 1)) == '1');

% Number of states:
N = size(X,1);

% Compute number of constraints:
m = 0;

for i = 1:K
    m = m + 2*n*nchoosek(n-1, i-1);
end

% Compute different constraints:
Os = zeros(N,n,m);
count = 1;

% Loop over variables:
for i = 1:n
    
    % Indices of states with variable i either 1 or 0:
    inds_up = (X(:,i) == 1);
    inds_down = (X(:,i) == 0);
    
    % For K = 1:
    Os(inds_up, i, count) = 1;
    count = count + 1;
    Os(inds_down, i, count) = 1;
    count = count + 1;
    
    % Loop over different orders of marginals for K > 1:
    for j = 2:K
        
        groups = nchoosek([1:(i-1),(i+1):n], j-1);
        
        % Loop over different groups of nodes to marginalize over:
        for l = 1:size(groups,1)
            
            group = groups(l,:);
            
            % Constrain correlation between group of nodes:
            Os(logical(inds_up.*prod(X(:,group), 2)), i, count) = 1;
            count = count + 1;
            Os(logical(inds_down.*prod(X(:,group), 2)), i, count) = 1;
            count = count + 1;
            
        end
    end
end
        
        