% Script to estimate local irreversibilities and compute minimum
% irreversibilities while correcting for finite-data effects. We loop over
% the different order of interactions from 1 to N (the number of cells in a
% group).

% Load state transitions:
data_struct = load('salamander_retina_transitions.mat');

% Transitions for each repeat of each stimulus and each cell group:
transitions = data_struct.transitions;

% Number of repeats for each stimulus:
num_repeats = data_struct.num_repeats;

% Lengths of different repeats of different stimuli (in ms):
stimuli_lengths = data_struct.stimuli_lengths;

% Number of stimuli and samples:
num_stimuli = data_struct.num_stimuli;
num_samples = data_struct.num_cellSamples;

% Size of cell groups:
n = data_struct.n;

%% Compute fractions of transitions for subsampled data:

% Data fractions to consider:
fracs = .9:-.1:.5;
num_fracs = length(fracs);

% Number of data subsamples for each data fraction:
num_dataSamples = 100;

% Initial step size for minimum irreversibility algorithm:
step_size_init = .5;

% Number of iterations of minimum irreversibility algorithm:
num_steps = 10^2;

% Transpose indices:
inds_trans = transpose_indices_multipartite(n);

% Observables for different constraints:
Os = cell(n-1,1);
for i = 1:(n-1)
    Os{i} = binary_constraints_multipartite(n, i);
end

% List of cell samples:
cell_samples = data_struct.cell_samples;

% Record number of transitions:
num_trans = zeros(num_stimuli, num_samples);

% Record irreversibilities and inifinite-data estimates:
irreversibilities = zeros(num_stimuli, num_samples, num_fracs + 1, num_dataSamples);
irreversibilities_inf = zeros(num_stimuli, num_samples, num_dataSamples);

% Record minimum irreversibilities of order k from 1 to N-1:
irreversibilities_min = zeros(num_stimuli, n-1, num_samples, num_fracs + 1, num_dataSamples);
irreversibilities_min_inf = zeros(num_stimuli, n-1, num_samples, num_dataSamples);

% Record proportions of minimum and interaction irreversibilities:
props_inf = zeros(num_stimuli, n-1, num_samples, num_dataSamples);
props_int_inf = zeros(num_stimuli, n-1, num_samples, num_dataSamples);

% Record changes in state probabilities:
change_stateProbs = zeros(num_stimuli, num_samples, 2^n, num_fracs + 1, num_dataSamples);
change_stateProbs_inf = zeros(num_stimuli, num_samples, 2^n, num_dataSamples);

% Loop over different stimuli:
for i = 1:num_stimuli
    
    % Loop over different group samples:
    for j = 1:num_samples
        
        tic
        
        % Copmute full transition matrix:
        T = zeros(2^n, n+1);
        T_shuffle = zeros(2^n, n+1);
        
        for k = 1:num_repeats(i)
            T = T + full(transitions{i,j}{k});
        end
        
        % Number of transitions:
        num_trans(i,j) = sum(T(:));
        
        % List all transitions:
        trans = [];
        
        for k = 1:2^n
            for l = 1:(n+1)
                trans = [trans; repmat([k, l], T(k,l), 1)];
            end
        end
        
        % Compute irreversibility with full data (using pseudocount method):
        P = (T(:,1:n) + 1)/sum(T(:) + 1);
        irreversibilities(i,j,1,:) = sum(P(:).*log2(P(:)./P(inds_trans)));
        
        % Compute independent irreversibility with full data:
        p_1 = reshape(sum(Os{1}.*repmat(P, 1, 1, 2*n), [1 2]), 2, n);
        irreversibilities_min(i,1,j,1,:) = sum(p_1.*log2(p_1./flipud(p_1)), [1 2]);
        
        % Compute change in state probabilities:
        P_full = convert_transitions_multipartite((T + 1)/sum(T(:) + 1));
        change_stateProbs(i,j,:,1,:) = repmat(reshape(sum(P_full,1) - sum(P_full,2)', 1, 1, 2^n), 1, 1, 1, 1, num_dataSamples);
        
        % Keep track of transitions:
        Ps = zeros(2^n, n, num_fracs, num_dataSamples);
        
        % Loop over data subsamples:
        for k = 1:num_dataSamples
        
            % Initialize list of transitions:
            trans_temp = trans;
            
            % Loop over fractions of data:
            for l = 1:num_fracs
                
                % Sample transitions:
                inds = randsample(size(trans_temp,1), round(fracs(l)*num_trans(i,j)), false);
                trans_temp = trans_temp(inds,:);
                [trans_unique, ~, ic] = unique(trans_temp, 'rows');
                T_temp = zeros(2^n, n + 1);
                for m = 1:size(trans_unique, 1)
                    T_temp(trans_unique(m,1), trans_unique(m,2)) = sum(ic == m);
                end
                
                % Compute irreversibility (using pseudocounts):
                P_temp = (T_temp(:,1:n) + 1)/sum(T_temp(:) + 1);
                irreversibilities(i,j,l+1,k) = sum(P_temp(:).*log2(P_temp(:)./P_temp(inds_trans)));
                
                Ps(:,:,l,k) = P_temp;
                
                % Compute independent irreversibility:
                p_1 = reshape(sum(Os{1}.*repmat(P_temp, 1, 1, 2*n), [1 2]), 2, n);
                irreversibilities_min(i,1,j,l+1,k) = sum(p_1.*log2(p_1./flipud(p_1)), [1 2]);
                
                % Compute change in state probabilities:
                P_full = convert_transitions_multipartite((T_temp + 1)/sum(T_temp(:) + 1));
                change_stateProbs(i,j,:,l+1,k) = reshape(sum(P_full,1) - sum(P_full,2)', 1, 1, 2^n);
                  
            end
            
            % Extrapolate irreversibility to infinite data with linear fit:
            fit = polyfit([1, fracs].^(-1), reshape(irreversibilities(i,j,:,k), 1, num_fracs + 1), 1);
            irreversibilities_inf(i,j,k) = fit(2);
            
            % Extrapolate independent irreversibility to infinite data:
            fit_ind = polyfit([1, fracs].^(-1), reshape(irreversibilities_min(i,1,j,:,k), 1, num_fracs + 1), 1);
            irreversibilities_min_inf(i,1,j,k) = fit_ind(2);
            fit_prop_ind = polyfit([1, fracs].^(-1), reshape(irreversibilities_min(i,1,j,:,k), 1, num_fracs + 1)./...
                reshape(irreversibilities(i,j,:,k), 1, num_fracs + 1), 1);
            props_inf(i,1,j,k) = fit_prop_ind(2);
            
            % Extrapolate to infinite data for changes in state probabilities:
            for l = 1:(2^n)
                fit_dP = polyfit([1, fracs].^(-1), reshape(change_stateProbs(i,j,l,:,k), 1, num_fracs + 1), 1);
                change_stateProbs_inf(i,j,l,k) = fit_dP(2);
            end
            
        end
        
        % If irreversibility is significant, then compute minimum irreversibilities of different orders:
        if mean(irreversibilities_inf(i,j,:)) > 2*std(irreversibilities_inf(i,j,:))
            
            % Compute minimum irreversibility for full data:
            P_temp = P;
            
            for k = (n-1):-1:2
                
                [P_temp, S_temp, ~] = min_irreversibility_multipartite(Os{k}, P_temp, num_steps, step_size_init);
                irreversibilities_min(i,k,j,1,:) = S_temp;
                
            end
            
            % Loop over data subsamples:
            for k = 1:num_dataSamples
                
                % Loop over fractions of data:
                for l = 1:num_fracs
                    
                    % Compute minimum irreversibility for full data:
                    P_temp = Ps(:,:,l,k);
                    
                    for m = (n-1):-1:2
                        
                        [P_temp, S_temp, ~] = min_irreversibility_multipartite(Os{m}, P_temp, num_steps, step_size_init);
                        irreversibilities_min(i,m,j,l+1,k) = S_temp;
                        
                    end
                    
                end
                
                % Extrapolate minimum irreversibilities to infinite data:
                for l = 2:(n-1)
                    
                    % Fit minimum irreversibility versus inverse data fraction:
                    fit_irreversibility = polyfit([1, fracs].^(-1), reshape(irreversibilities_min(i,l,j,:,k), 1, num_fracs + 1), 1);
                    irreversibilities_min_inf(i,l,j,k) = fit_irreversibility(2);
                    
                    % Fit minimum irreversibility proportion versus inverse data fraction:
                    fit_prop = polyfit([1, fracs].^(-1), reshape(irreversibilities_min(i,l,j,:,k), 1, num_fracs + 1)./...
                        reshape(irreversibilities(i,j,:,k), 1, num_fracs + 1), 1);
                    props_inf(i,l,j,k) = fit_prop(2);
                    
                    % Fit interaction irreversibility proportion versus inverse data fraction:
                    fit_prop_int = polyfit([1, fracs].^(-1), reshape(irreversibilities_min(i,l,j,:,k) - irreversibilities_min(i,l-1,j,:,k), 1, num_fracs + 1)./...
                        reshape(irreversibilities(i,j,:,k), 1, num_fracs + 1), 1);
                    props_int_inf(i,l,j,k) = fit_prop_int(2);
                    
                end
            end
        end
        
        % Print some things:
        i
        j
        toc
        
    end
end

% First-order interaction irreversibilities are same as first-order minimum irreversibilities:
props_int_inf(:,1,:,:) = props_inf(:,1,:,:);

% Compute significance of different cell groups during different stimuli:
sig = (abs(mean(irreversibilities_inf, 3)) > 2*std(irreversibilities_inf, [], 3)).*sign(mean(irreversibilities_inf, 3));
sig_min = (abs(mean(irreversibilities_min_inf, 4)) > 2*std(irreversibilities_min_inf, [], 4)).*sign(mean(irreversibilities_min_inf, 4));
sig_change_stateProbs = (abs(mean(change_stateProbs_inf, 4)) > 2*std(change_stateProbs_inf, [], 4)).*sign(mean(change_stateProbs_inf, 4));

% Indices of significant cell groups:
inds_sig = cell(num_stimuli,1);
inds_sig_min = cell(num_stimuli, n-1);

for i = 1:num_stimuli
    
    inds_sig{i} = find(sig(i,:) > 0);
    
    for j = 1:(n-1)
        inds_sig_min{i,j} = find(sig_min(i,j,:) > 0);
    end
end

% Save results:
save('salamander_retina_irreversibility', 'num_stimuli', 'num_repeats',...
    'num_dataSamples', 'fracs', 'stimuli_lengths', 'cell_samples', 'num_trans',...
    'irreversibilities', 'irreversibilities_min', 'irreversibilities_inf',...
    'irreversibilities_min_inf', 'props_inf', 'props_int_inf', 'sig', 'sig_min',...
    'inds_sig', 'inds_sig_min');


