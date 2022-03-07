% Script to compute the transitions between binary states for
% randomly-selected groups of cells. For each cell group, we record the
% state transitions for each repeat of each stimulus.

% Load spiking:
data_struct = load('Data_processed.mat');

% Pick stimuli:
stimuli = 1:3;

% Cell IDs for spikes responding to different stimuli:
cell_IDs = data_struct.cell_IDs(stimuli);

% Spike times in different stimuli:
spike_times = data_struct.spike_times(stimuli);

% Lengths of each stimulus repeat (in ms):
stimuli_lengths = data_struct.stimuli_lengths(stimuli);

% Numbers of cells, stimuli, and repeats:
num_cells = data_struct.num_cells;
num_stimuli = length(stimuli);
num_repeats = data_struct.num_repeats(stimuli);

% Choose group size:
n = 5;

%% Compute transitions between states for randomly-generated cell groups of a given size:

% Number of random cell groups:
num_cellSamples = 100;

% Window size (in ms):
dt = 20;

% Things to record:
cell_samples = cell(1, num_cellSamples);
transitions = cell(num_stimuli, num_cellSamples);

% Loop over different group samples:
for i = 1:num_cellSamples
    
    % Pick random group of cells:
    cells = datasample(1:num_cells, n, 'Replace', false);
        
    % Record group of cells:
    cell_samples{i} = cells;
    
    % Loop over different stimuli:
    for j = 1:num_stimuli
        
        % Record transitions:
        transitions_temp = cell(num_repeats(j), 1);
        
        % Loop over different stimulus repeats:
        for k = 1:num_repeats(j)
            
            % Compute transitions:
            T = transitions_slidingWindow_variableLengths(spike_times{j}(k),...
                cell_IDs{j}(k), cells, dt, stimuli_lengths{j}(k));
            transitions_temp{k} = sparse(T);
            
        end
        
        % Record things:
        transitions{j,i} = transitions_temp;
        
    end
end

% Save data:
save('salamander_retina_transitions', 'n', 'stimuli', 'num_stimuli',...
    'num_cellSamples', 'num_repeats', 'stimuli_lengths', 'cell_samples', 'transitions');
                