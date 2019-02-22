function [SSIMS_coordinates] = getSSIMSprojection(sorted_timestamps, ev, starttime, stduration, q, tSNE_transform, basespiketrains)

% [SSIMS_coordinates] = getSSIMSprojection(sorted_timestamps, ev, starttime, stduration, q, tSNE_transform, basespiketrains)
%
% Usage:
% Once you have generated a Spike train SIMilarity Space (use getSSIMS function)
% use this function to project new data onto it.
%
% INPUTS:
% sorted_timestamps: cell array of neuron spiking stimestamps
% ev: events to align to
% start_times: offsets relative to event, vector of times at which projections will be calculated
% stduration: length of windows to consider
% q: parameter for spike train metric algorithm (1/q ~temporal accuracy)
% tSNE_transform: transform to apply to distance matrices (from getSSIMS)
% basespiketrains: spike trains that make up the base of the space (from getSSIMS)
%
% OUTPUT:
% SSIMS_coordinates: n-dim coordinates of the points associated with each spike train, where
% 	n = size(tSNE_transform, 2);
%

% @author Carlos Vargas-Irwin
% Copyright (c) Carlos Vargas-Irwin, Brown University. All rights reserved.

warning('SSIMSToolBox:getSSIMSprojection:noCompiledFunction', 'Using plain MATLAB implementation of the SSIMS algorithm. This calculation will be very slow. For speed gains of 1-2 orders of magnitude, consider compilation (see doc/INSTALL.md).');


%%%%%%%%%%%%  Calculate Neural Trajectories  %%%%%%%%%%%%%%%%%%%

%multiWaitbar('CloseAll');
multiWaitbar('Timesteps', 0);
multiWaitbar('Unit', 'Color', [.5, .2, 0]);
multiWaitbar('Timesteps', 'color', [.5, 0,0]);

dim = size(tSNE_transform, 2);
base_l = size(tSNE_transform, 1);
n_n = numel(sorted_timestamps);
n_base_ev = base_l / n_n;
n_ev = numel(ev);

SSIMS_coordinates = zeros (length(ev),dim,length(starttime));
tic
for st = 1:length(starttime)

    multiWaitbar('Timesteps', 'Value', st/length(starttime));
    multiWaitbar('Unit', 0);

    spikes = spikeTrainEventWindowSpikes(sorted_timestamps, ev + starttime(st), 0, stduration);
    dmat = zeros(n_ev, base_l);

    for i_n = 1:length(sorted_timestamps)
        multiWaitbar('Unit', 'Value', i_n/length(sorted_timestamps));
        dmat(:, (1:n_base_ev) + n_base_ev * (i_n - 1) )  = spikeTrainDistMat(spikes{i_n}, basespiketrains{i_n}, q);
    end
    SSIMS_coordinates(:, :, st) = dmat * tSNE_transform;
end

telapsed = toc;
display (sprintf('Calculated the projections for %i units and %i timesteps in %4.1f seconds.\n', n_n, length(starttime), telapsed ))

multiWaitbar('Timesteps', 'Close');
multiWaitbar('Unit', 'Close');

