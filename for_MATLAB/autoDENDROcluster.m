function [s   cindex     nclus   CI] = autoDENDROcluster(crange,datamat)
%%  SPIKE TRAIN SIMILARITY SPACE TOOLBOX% 
% [s cindex centroids nclus ] = autokmeanscluster(crange,datamat)
% helper function for SSIMSETS
%
% @author Carlos Vargas-Irwin
% Copyright (c) Carlos Vargas-Irwin, Brown University. All rightsreserved. Resistance is futile.


% calculate silhouette value for number of clusters specified in crange
 L = linkage(datamat);

for c = 2:max(crange)
     
    cix{c} = cluster(L,'maxclust',c);  
 
[sv(:,c) ] = silhouette(datamat,cix{c});

[s(c)] = mean(sv(:,c));
[  CI(c,:)  ] =  prctile(sv(:,c), [5 95]);
 

% [ mu STD ci(c,:) ] = normfit( sv(:,c) );
            
end

% cluster according to mas silhouette value

[v nclus] =  max(s);
%nclus=4; 
cindex = cix{nclus};
 

