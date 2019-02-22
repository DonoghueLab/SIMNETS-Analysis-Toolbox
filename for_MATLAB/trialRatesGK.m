 function [frmat] = trialRatesGK(sorted_timestamps,event,edges,gkbins)
% [frmat] = trialRatesGK(sorted_timestamps,event,edges,gkbins)
%
% INPUTS:
%  sorted_timestamps: cell array of neuron spiking stimestamps
%  event: events to align to
%  edges: bin edges relative to event (e.g. [0:0.001:1.5]
%  gkbins: number of bins to use for gaussian kernel 
%
% OUTPUT:
% FRmat: matrix of firing rates smoothed with gaussian kernel (size =  #events x #bins x #neurons)


%get spike trains
spk = spikeTrainEventWindowSpikes(sorted_timestamps, event, 0, edges(end));

gaussk = gausswin(gkbins);

% convolve with gaussian kernel to get continuous firing rates
frmat = zeros(numel(event),numel(edges)-1,numel(sorted_timestamps));
for n = 1:numel(sorted_timestamps)
    for e = 1:numel(event)
        s1 = spk{n}{e};
        c = histc(s1,edges);
        c = c(1:end-1);
        cs = conv(c, gaussk,'same');
        if(isempty(cs))
            cs = zeros(1,numel(edges)-1) ;
        end
        frmat(e,:,n) = cs;
    end
end
