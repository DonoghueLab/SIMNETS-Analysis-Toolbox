function [distm] = spikeTrainDistMat(spikes1, spikes2, q)
% 	SPIKETRAINDISTMAT   Compiles distance matrix for spike trains
% 		[DISTM] = SPIKETRAINDISTMAT(SPIKES1, SPIKES2, Q)
%
% 	Given two cells of spike trains and a cost parameter, this function compiles the distance matrix
%	according to the Victor and Purpura metric.
%
% 	Created by Jonas B Zimmermann on 2013-09-13.
% 	Copyright (c) 2013 Donoghue Lab, Neuroscience Department, Brown University.
% 	All rights reserved.


distm = zeros(numel(spikes1), numel(spikes2));
for c = 1:numel(spikes1)
	for r = 1:numel(spikes2)
		distm(c, r) = vpSpikeTimeDist(spikes1{c}, spikes2{r}, q);
	end
end

end %  function
