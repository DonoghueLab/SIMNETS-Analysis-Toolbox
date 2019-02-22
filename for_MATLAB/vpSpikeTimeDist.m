function [d dmat] = vpSpikeTimeDist(s1,s2,q)
% Spike train metric as described by Victor and Purpura
% @author Carlos Vargas-Irwin
% Copyright (c) Carlos Vargas-Irwin, Brown University. All rights reserved.

dmat = zeros(length(s1)+1, length(s2)+1);

dmat(:,1) = 0:length(s1);
dmat(1,:) = 0:length(s2);

for r = 2:length(s1)+1
    for c = 2:length(s2)+1
        d = abs( s1(r-1) - s2(c-1)  );
        dmat(r,c) =  min( [ dmat(r,c-1)+1    dmat(r-1,c)+1    dmat(r-1,c-1)+(q*d)   ] );
    end
end

d = dmat(length(s1)+1, length(s2)+1);

