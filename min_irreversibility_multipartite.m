function [P, S, S0, count] = min_irreversibility_multipartite(Os, P0, num_steps, step_size_init)
% Inputs: 2^n x n x m array of constraints Os, where Os(:,:,i) is the ith
% constraint and n is the number of variables in the system; 2^n x n matrix
% of initial transition probabilities P0, where P0(i,j) represents variable
% j changing in state i (we note that sum(P0(:)) <= 0 and P0 >= 0); number
% of iterations of algorithm num_steps; and initial step size
% step_size_init.
%
% Outputs: 2^n x n transition probability matrix P that minimizes
% irreversibility while matching the constraints; final irreversibility S;
% initial irreversibility S0; and number of iterations of minimization
% algorithm.
%
% NOTE: This is an adaptation of the Frank-Wolfe algorithm.

% Threshhold for change in model:
dP_threshold = 10^(-10);

% Number of variables:
n = size(Os,2);

% Number of constraints:
m = size(Os,3);

% Options for linear program solver:
options = optimoptions('linprog', 'Display', 'none', 'ConstraintTolerance', 10^(-9));

% Get transpose indices:
inds_trans = transpose_indices_multipartite(n);

% Zero out transitions with zero reverse-probability:
P0(P0(inds_trans) == 0) = 0;

% Only keep non-zero transition probabilities:
inds_pos = find(P0(:) > 0);
P0 = P0(inds_pos);
[~, inds_trans] = sort(inds_trans(inds_pos));
N = length(inds_pos);

% Compute initial irreversibility:
S0 = sum(P0.*log2(P0./P0(inds_trans)));

% Turn each constraint into an array:
Os_lin = zeros(N, m);
for i = 1:m
    Os_temp = Os(:,:,i);
    Os_lin(:,i) = Os_temp(inds_pos);
end

% Remove constraints with zero entries:
Os_lin = Os_lin(:, sum(abs(Os_lin),1) > 0);
m = size(Os_lin,2);

% Compute average operators:
Os_bar = zeros(m,1);
for i = 1:m
    Os_bar(i) = sum(P0.*Os_lin(:,i));
end

% If P0 is symmetric then return:
if sum(abs(P0 - P0(inds_trans))) < dP_threshold
    P = zeros(2^n, n);
    P(inds_pos) = P0;
    S = S0;
    count = 0;
    return;
end

% Initialize variables to return:
P = P0;
S = S0;

% Loop until the model changes less than the threshold:
P_new = P0;

for i = 1:num_steps

    % Record old transition probabilities:
    P_old = P_new;
    
    % Compute derivatives of irreversibility:
    dS_dP = -log(P_old(inds_trans)./P_old) - P_old(inds_trans)./P_old + 1;
    
    % Solve linear program for local optimization:
    P_temp = linprog(dS_dP, ones(1,N), 1, Os_lin', Os_bar, zeros(N,1), ones(N,1), options);
    
    % Return if linear program solver fails:
    if min(P_temp) < 0
        break;
    end
    
    % Choose step size:
    step_size = step_size_init/i;
    
    % New transition probabilities:
    P_new = P_old + step_size*(P_temp - P_old);
    
    % Compute change in model:
    dP = sum(abs(P_new - P_old));
    
    % Keep track of minimum irreversibility:
    S_new = sum(P_new.*log2(P_new./P_new(inds_trans)));
    if S_new < S
        P = P_new;
        S = S_new;
    end
    
    % Break from loop if change in model is smaller than threshold:
    if dP < dP_threshold
        break;
    end
    
end

% Finalize variables to return:
P_temp = P;
P = zeros(2^n, n);
P(inds_pos) = P_temp;
count = i;



