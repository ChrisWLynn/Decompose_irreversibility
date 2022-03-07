function inds = transpose_indices_multipartite(n)
% Input: Number of binary variables n
%
% Output: Consider a 2^n x n transition matrix P, where P(i,j) indicates
% that variable j flips in state i. We return the linear indices inds, such
% that reshape(P(inds), 2^n, n) is the time-reverse of P

% Matrix of states:
X = double(dec2bin(0:(2^n - 1)) == '1');

% Number of states:
N = size(X,1);

% Initialize indices:
inds = zeros(N*n, 1);

% Loop over states:
for i = 1:N
    
    x = X(i,:);
    
    % Loop over variables:
    for j = 1:n
        
        y = x;
        y(j) = double(~x(j));
        ind_y = bin2dec(strjoin(string(y),'')) + 1;
        
        % Add transpose index:
        inds(sub2ind([N, n], i, j)) = sub2ind([N, n], ind_y, j);
        
    end
    
end

