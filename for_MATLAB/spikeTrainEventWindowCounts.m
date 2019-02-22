function [counts] = spikeTrainEventWindowCounts(spkT, events, offsets, lengths, varargin)
% 	SPIKETRAINEVENTWINDOWCOUNTS   returns spike counts for spike trains, at events around offsets in windows of length lengths
% 		[COUNTS] = SPIKETRAINEVENTWINDOWCOUNTS(SPKT, EVENTS, OFFSETS, LENGTHS[, NORMALIZE])
%
% 	Arguments:
% 		spkT is a single spike train (vector of spike times [double]), or a cell array of spike trains
% 		events is a vector of event times (double) or a cell array of event time vectors
% 		offsets is a single offset time (relative to event times) for spike train window, or it is a vector with numel(events) elements (only if events is cell)
% 		lengths is a single window length (relative to offsets) for spike train window, or it is a vector with numel(events) elements (only if events is cell)
%		normalize: if true, spike counts will be normalized by window length {0, [1]}
%	Returns:
%		count: cell array; numel(spkT) x numel(events) (assuming spkT and events being cell arrays)
% 	Created by Jonas B Zimmermann on 2015-10-28.
% 	Copyright (c) 2015 Donoghue Lab, Neuroscience Department, Brown University.
% 	All rights reserved.
% 	@author  $Author: Jonas B Zimmermann$

if ~iscell(spkT)
	spkT = {spkT};
end
n_spkt = numel(spkT);

if ~iscell(events)
	events = {events};
end
n_events = numel(events);

if numel(offsets) == 1
	offsets = repmat(offsets, n_events, 1);
else
	if numel(offsets) ~= n_events
		error('SSIMS:spikeTrainEventWindowCounts:wrongNOffsets', 'Argument ''offsets'' has to have the same number of elements as ''events'' or it has to be scalar.');
	end
	offsets = reshape(offsets, n_events, []);
end

if numel(lengths) == 1
	lengths = repmat(lengths, n_events, 1);
else
	if numel(lengths) ~= n_events
		error('SSIMS:spikeTrainEventWindowCounts:wrongNOffsets', 'Argument ''lengths'' has to have the same number of elements as ''events'' or it has to be scalar.');
	end
	lengths = reshape(lengths, n_events, []);
end

if nargin < 5
	normalize = 1;
else
	normalize = varargin{1};
end


counts = cell(n_spkt, n_events);

for i_s = 1:n_spkt
	for i_e = 1:n_events
		if normalize
			n_factor = 1.0 / lengths(i_e);
		else
			n_factor = 1.0;
		end
		n_this_event = numel(events{i_e});
		counts{i_s, i_e} = zeros(n_this_event, 1);
		for i_t = 1:n_this_event
			counts{i_s, i_e}(i_t) = sum((events{i_e}(i_t) + offsets(i_e) <= spkT{i_s}) & (spkT{i_s} < events{i_e}(i_t) + offsets(i_e) + lengths(i_e))) * n_factor;
		end
	end
end

end %  function
