function T = transitions_slidingWindow_variableLengths(spike_times, cell_IDs, cells, dt, trial_lengths)
% Input: cell array of spike times (in ms) for different stimulus repeats
% (where spike_times{i} gives the array of spike times in repeat i), cell
% array of cell IDs in order of spiking (where cell_IDs{i} gives the array
% of cells that spiked in repeat i), the subset of cells that we want to
% study cells, the width of the time window dt (in ms), and array
% trial_lengths where trial_lengths(i) is the length of trial i (in ms).
%
% Output: Calculate state transitions in the following way. We move a
% sliding window of width dt and record the transitions between binary
% states. For times with simultaneous spikes, we pick a random order. We
% return the 2^n x n matrix of transition counts T, where n is the number
% of cells. For example, T(i,j) is the number of times cell j flipped when
% in state i.

% Number of trials:
num_trials = length(spike_times);

% Number of cells:
n = length(cells);

% Matrix of states:
num_states = 2^n;

% Create indices for cells:
cell_inds(cells) = 1:n;

% Initialize transition matrix:
T = zeros(num_states, n+1);

% Loop over trials:
for i = 1:num_trials
    
%     tic
    
    % Spiking data for this trial:
    times_trial = spike_times{i};
    cells_trial = cell_IDs{i};
    
    % Restrict to the desired set of cells:
    inds = ismember(cells_trial, cells);
    times_trial = times_trial(inds);
    cells_trial = cell_inds(cells_trial(inds));
    
    % Number of spikes in this trial:
    num_spikes = length(times_trial);
    
    % Add small noise to break simultaneous spikes:
    times_trial = times_trial + 10^(-6)*(rand(1, num_spikes) - .5);
    [~, order] = sort(times_trial);
    times_trial = times_trial(order);
    cells_trial = cells_trial(order);
    
    % Initialize first state in this trial:
    inds = find(times_trial < dt);
    times_temp = times_trial(inds);
    cells_temp = cells_trial(inds);
    
    state_old = zeros(n,1);
    state_old(cells_temp) = 1;
    ind = bin2dec(join(string(state_old),'')) + 1;
    
    % Loop over remaining spikes:
    for j = (length(inds) + 1):num_spikes
        
        time_up = times_trial(j);
        cell_up = cells_trial(j);
        
        % Number of spikes to remove before adding new spike:
        num_remove = sum(times_temp < time_up - dt);
        
        % Loop over spikes to remove:
        for k = 1:num_remove
            
            % Remove spike:
            times_temp = times_temp(2:end);
            cells_temp = cells_temp(2:end);
            
            state_new = zeros(n,1);
            state_new(cells_temp) = 1;
            
            % Node that flipped:
            node = find(abs(state_new - state_old));
            if isempty(node)
                node = n + 1;
            end
            
            % Increment transition count:
            T(ind, node) = T(ind, node) + 1;
            
            % Update system state:
            state_old = state_new;
            ind = bin2dec(join(string(state_old),'')) + 1;
            
        end
        
        % Add new spike:
        times_temp = [times_temp, time_up];
        cells_temp = [cells_temp, cell_up];
        
        state_new = zeros(n,1);
        state_new(cells_temp) = 1;
        
        % Node that flipped:
        node = find(abs(state_new - state_old));
        if isempty(node)
            node = n + 1;
        end
        
        % Increment transition count:
        T(ind, node) = T(ind, node) + 1;
        
        % Update system state:
        state_old = state_new;
        ind = bin2dec(join(string(state_old),'')) + 1;
        
    end
    
    % Remove remaining spikes unitl we hit the end of the trial:
    num_remove = sum(times_temp < trial_lengths(i) - dt);
    
    % Loop over spikes to remove:
    for j = 1:num_remove
        
        times_temp = times_temp(2:end);
        cells_temp = cells_temp(2:end);
        
        state_new = zeros(n,1);
        state_new(cells_temp) = 1;
        
        % Node that flipped:
        node = find(abs(state_new - state_old));
        if isempty(node)
            node = n + 1;
        end
        
        % Increment transition count:
        T(ind, node) = T(ind, node) + 1;
        
        % Update system state:
        state_old = state_new;
        ind = bin2dec(join(string(state_old),'')) + 1;
        
    end
    
end
    
