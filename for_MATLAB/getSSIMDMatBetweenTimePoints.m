
function [dmat, basespiketrains] = getSSIMDMatBetweenTimePoints(sorted_timestamps, ev_list, stduration, q)
% [dmat, basespiketrains] = getSSIMDMatBetweenTimePoints(sorted_timestamps, ev_list, stduration, q)
%  SPIKE TRAIN SIMILARITY SPACE TOOLBOX
%
% OVERVIEW:
%  This function generates a spike train distance matrix.
%
%
% INPUTS: (all times in seconds)
% sorted_timestamps = cell array with the sorted timestamps for each neuron (length = n)
% ev_list = timestamps for events to align spike trains to (length = m)
% stduration = duration of spike trains to be analyzed
% q = from Victor & Purpura's algorithm, 1/q specifies temporal precision
%
% OUTPUTS:
% dmat = distance matrix of size m×(n m);
%   for each neuron an m×m matrix will be calculated, and those are concatenated column-wise
% basespiketrains = cell array of spike trains used as the basis for the space (relative to specified events)
%
% @author Carlos Vargas-Irwin
% Copyright (c) Carlos Vargas-Irwin, Brown University. All rights reserved.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basespiketrains = spikeTrainEventWindowSpikes(sorted_timestamps, ev_list, 0, stduration);

n_ev = numel(ev_list);
dmat = zeros(n_ev, numel(sorted_timestamps) * n_ev);

for i_n = 1:numel(sorted_timestamps)
    
    dmat(:, (1:n_ev) + n_ev * (i_n - 1) ) = spikeTrainDistMat(basespiketrains{i_n}, basespiketrains{i_n}, q);
    
end
