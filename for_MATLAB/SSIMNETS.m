function [ SNETS  ] = SSIMNETS(sorted_timestamps, event, stduration, q, perp, varargin)
%  SPIKE TRAIN SIMILARITY SPACE TOOLBOX% 
%  OVERVIEW:
%  This function identifies 'functional sub-ensembles' (sub-nets)
%  This is done by clustering neurons according to the correlation between
%  single unit SSIMS relational maps
%  NOTE: the firing patterns are not compared directly! Only in terms of
%  their SSIMS representation. So: neurons that rank the same trials as
%  being similar will be clustered together (even if their firing patterns
%  are different.
%  Neurons are clustered using K-means, the number of clusters is selected 
%  using silhouette values 
%  (high values ~ higer between/within cluster distances)
%
% INPUTS: (all times in seconds)
% sorted_timestamps = cell array with the sorted timestamps for each neuron (length = n)
% event = timestamps for events to align spike trains to (length = m)
% stduration = duration of spike trains to be analyzed
% q = from Victor & Purpura's algorithm, 1/q specifies temporal precision
% neurondim  = # of dimensions used to cluster neurons
% perp = perplexity used to cluster neurons
% crange = # of clusters to check using k means
% displayopt = the last (optional) argument will plot results if set to 1
%
% Code requires SSIMS toolbox + helper functions:
%   - SSIMSNeuronRelCorr
%   - autokmeanscluster
%
% OUTPUTS: 
% SNETS.NSPACE: each point represents a neuron
% SNETS.clusterindex: index denoting membership in a cluster
% SNETS.centroids: cluster centroids
% SNETS.silhouette: silhouette values for different # of clusters identified using k-means (Higher = more cluster separation)
% SNETS.numclus: # of clusters identified (NOTE: this version will only detect > 1 cluster!!!)
% SNETS.NNcorr: Correlations between single unit SSIMS relational maps ( neuron x neuron )
% SNETS.distmat: full distance matrix ( neuron x events )
%
% 
% Questions? Contact Carlos_vargas_Irwin@brown.edu or Jacqueline Hynes@Brown.edu
% @author Jacqueline Hynes
% Copyright (c) Jacqueline Hynes, Brown University. All rights reserved.
%% Defaults:

p = inputParser;
                         
optionalsDefault = {3, 3, 0, 10, 99,'pearson','tsne'};                    % Placeholder for default values

checknum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
checkOn = @(x) isnumeric(x) && isscalar(x);

addParameter(p,'displaydim',optionalsDefault{1},checknum)
addParameter(p,'clusterdim',optionalsDefault{2},checknum)
addParameter(p,'iterations',optionalsDefault{3},checkOn)
addParameter(p,'crange',optionalsDefault{4},checknum)
addParameter(p,'silhCI',optionalsDefault{5},checknum)
addParameter(p,'corrType',optionalsDefault{6}, @ischar)
addParameter(p,'drTechnique',optionalsDefault{7}, @ischar)

parse(p,varargin{:})
varOption =p.Results;
 
%% STEPS 1-2: Calculate the SSIM matrices for each neuron, the pariwise correlations between them
 
 % Generate the NxN neuron similarity matrix with or without shuffle-based
 % statistical test for sub-net (cluster) number validation  

 [ncmat, dmatnb, basespiketrains] = SSIMSNeuronRelCorr(sorted_timestamps, event, stduration, q, varOption.corrType); 
 

%% STEPS 3: Dimensionality Reduction
% Use a dimensionality reduction method to project the matrix of
% neuron-neuron correlations into
% a lower dimensional space for clustering
  
switch varOption.drTechnique
        
    case 'pca'
        
    pcY = pca(ncmat);
    NSPACEcluster = ncmat * pc(:,1:varOption.clusterdim);

    case 'mds'

    mdsY = cmdscale(1-ncmat);   
    NSPACEcluster = neuralSpace(:,1:varOption.clusterdim);

    case 'tsne'

    NSPACEcluster = runtSNE(ncmat,varOption.clusterdim,perp); 
end

 [LPC] = pca(NSPACEcluster);
 NSPACE = NSPACEcluster*LPC(:,[1:varOption.displaydim]);
 
%% Step 4: Cluster detection and visulization 
% using K-means, selecting the number of clusters

 [silh, cindex,~  ,nclus] = autokmeanscluster(p.Results.crange,NSPACEcluster);  
 
 % Calculate significance of Silhouette Values
 iterations = varOption.clusterdim;
 
 [  ncmatShift, silhDistMU, silhDistCI] = SSIMSNeuronRelCorr_shuffleTest( dmatnb, iterations,perp, varOption.clusterdim, varOption.crange,varOption.silhCI,varOption.corrType,varOption.drTechnique );

 [ vals vIndex] = sort(silh, 'descend'); 
 
  if vals(1)<silhDistCI(vIndex(1))  % if the Max Avg. Silhouette value is not significant.. 

     [ cindex   ] = selectSignificanceClusterNumber(silh, silhDistCI, p.Results.crange,NSPACEcluster);
     nclus = unique(cindex); 
      
  end 


%% STORE DATA

    SNETS.NSmap = NSPACE;
    SNETS.clusterindex = cindex;
    SNETS.silh = silh;
    SNETS.shufMU = silhDistMU;
    SNETS.shufCI = silhDistCI;
    SNETS.numclus = nclus;
    SNETS.NNcorr = ncmat;
    SNETS.NNcorrS = ncmatShift;
    SNETS.distmat = dmatnb;
    SNETS.spiketrains = basespiketrains; 
    

end






