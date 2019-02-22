function [pcp] = getPCP(waves, V, npc)
%[pcp] = getPCP(waves, V, npc)
zm_wf = zeroMean(waves);

pcp = (V(:,1:npc)'*zm_wf);



