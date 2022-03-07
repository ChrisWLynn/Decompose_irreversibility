function P_full = convert_transitions_multipartite(P)
% Input: 2^n x n+1 multipartite transition matrix P, where n is the number
% of binary variables in the system and P(i,j) is the probability that 
% variable j changes in state i.
%
% Output: The full 2^n x 2^n transition matrix between states P_full

% Number of variables and states:
n = size(P,2) - 1;
N = size(P,1);

% Matrix of states:
X = double(dec2bin(0:(2^n - 1)) == '1')';

% Initialize full transition matrix:
P_full = zeros(N);

% Self-transition probabilities:
P_full(logical(eye(N))) = P(:,end);

% Loop over states:
for i = 1:N
    
    % Current state:
    x = X(:,i);
    
    % Loop over variables to flip:
    for j = 1:n
        
        % Next state:
        y = x;
        y(j) = double(~x(j));
        ind_y = bin2dec(strjoin(string(y),'')) + 1;
        
        % Fill in transition matrix:
        P_full(i, ind_y) = P(i,j);
        
    end
end

