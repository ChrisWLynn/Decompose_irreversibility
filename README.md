# Decompose_irreversibility
Code and data used to perform the analyses in the papers "Decomposing local irreversibility in interacting systems" and "Emergence of local irreversibility in complex interacting systems" by Christopher W. Lynn, Caroline M. Holmes, William Bialek, and David J. Schwab.

The processed neuronal spiking data is in "Data_processed.mat". The stimuli are in the following order: (1) Natural movie, (2) Brownian bar, and (3) repeated Brownian bar (the same trajectory is repeated multiple times). This data is published in O Marre, D Amodei, N Deshmukh, K Sadeghi, F Soo, TE Holy, and M J Berry II, “Mapping a complete neural population in the retina,” J. Neurosci. 32, 14859–14873 (2012) and S E Palmer, O Marre, M J Berry II, and W Bialek, “Predictive information in a sensory population,” Proc. Natl. Acad. Sci. 112, 6908–6913 (2015).

In order to estimate local irreversibilities and estimate kth-order minimum irreversibilities, one should first run "script_transitions.m". This script selects randomly-sampled groups of neurons and computes the transitions between binary states. One should then run "script_minimize_irreversibility.m". This script estimates the local irreversibility and kth-order minimum irreversibilities for each group of neurons and each stimulus while correcting for finite-data effects.

The above scripts use the following helper functions:

"transitions_slidingWindow_variableLengths.m" computes the transitions between binary states from neuronal spiking data.

"min_irreversibility_multipartite.m" minimizes irreversibility while preserving desired constraints on the system dynamics.

"binary_constraints_multipartite.m" gives the constraints for kth-order dynamics in binary multipartite systems.

"convert_transitions_multipartite.m" converts between different representations of binary multipartite transitions. 

"transpose_indices_multipartite.m" gives indices to transpose binary multipartite transitions, yielding reverse-time transitions.
