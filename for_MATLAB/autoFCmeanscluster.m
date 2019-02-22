function [s   cindex   centroids nclus   CI] = autoFCmeanscluster(crange,datamat)
%%  SPIKE TRAIN SIMILARITY SPACE TOOLBOX% 
% [s cindex centroids nclus ] = autokmeanscluster(crange,datamat)
% helper function for SSIMSETS
%
% @author Carlos Vargas-Irwin
% Copyright (c) Carlos Vargas-Irwin, Brown University. All rightsreserved. Resistance is futile.


% calculate silhouette value for number of clusters specified in crange
for c = 2:max(crange)
[cix{c} centc{c}] = kmeans(datamat, c,'Replicates',10);
[sv(:,c) ] = silhouette(datamat,cix{c});

[s(c)] = mean(sv(:,c));
[  CI(c,:)  ] =  prctile(sv(:,c), [5 95]);
 

% [ mu STD ci(c,:) ] = normfit( sv(:,c) );
            
end

% cluster according to mas silhouette value
[v nclus] =  max(s);
cindex = cix{nclus};
centroids = centc{nclus};

%% FCM 
% calculate silhouette value for number of clusters specified in crange

for c = 2:max(crange)
    
    [center{c} U obj_fcn] = fcm(datamat, c );
    
    maxU = max(U);
    clusterindex{c}(find(U(1,:) == maxU))=1;
    clusterindex{c}(find(U(2,:) == maxU))=2;
    [sv(:,c) ] = silhouette(datamat,clusterindex{c});

    [s(c)] = mean(sv(:,c));
    [  CI(c,:)  ] =  prctile(sv(:,c), [5 95]);

    % [ mu STD ci(c,:) ] = normfit( sv(:,c) );
            
end
[v nclus] =  max(s);

 cindex= clusterindex{nclus};
 centroids = center{nclus};