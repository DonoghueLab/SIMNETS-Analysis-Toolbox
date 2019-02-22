function [ ncmatShift, shufMU, shufCI] = SSIMSNeuronRelCorr_shuffleTest( dmatnb, it, perp, varargin );
%% SPIKE TRAIN SIMILARITY SPACE (SSIMS) AND SIMILARITY NETWORKS (SIMNETS) TOOLBOX% 
% OVERVIEW:
%
%
% INPUTS: 

%
%
% OUTPUTS: 
%
%
%
% EXAMPLE: [ ncmat, dmatnb, basespiketrains] = SSIMSNeuronRelCorr(sorted_timestamps, ev, 1, 10, 'pearson')
%
% Validated using MatLab R2016a,R2017a,R2018a, 
% See Papers for Details on SSIMS and SIMNETS Methods: 
%
% Jacqueline Hynes, David Brandman,  Jonas Zimmerman, John Donoghue, Carlos Vargas-Irwin (2018) 
% SIMNETS: a computationally efficient and scalable framework for identifying networks of functionally similar neurons
% 
% Vargas-Irwin, C. E., Brandman, D. M., Zimmermann, J. B., Donoghue, J. P., & Black, M. J
% "Spike Train SIMilarity Space (SSIMS): A Framework for Single Neuron and Ensemble Data Analysis (2014)."
% 
% Questions? Contact Carlos_vargas_Irwin@brown.edu or Jacqueline Hynes@Brown.edu
% @author Jacqueline Hynes
% Copyright (c) Jacqueline Hynes, Brown University. All rights reserved.

%% Optional arguments: 

optionalsDefault = { 10 , 10 , 99 , 'pearson', 'tsne'}; 
numInputs =  numel(varargin); % a and b are required

inputVar = 1;
if numInputs > 0

    while numInputs > 0
        if ~isempty(varargin{inputVar})
            optionalsDefault{inputVar} = varargin{inputVar};
        end
        inputVar = inputVar + 1;
        numInputs = numInputs - 1;
    end
end  
clusterdim = optionalsDefault{1}; 
crange = optionalsDefault{2};
silhCI = optionalsDefault{3};
corrType = optionalsDefault{4};
drTechnique = optionalsDefault{5};
 
%% Calculate pairwise distance matrix for each neuron
% Insert Explanation Here: Genertates N contatenated SxS Spike train Similarity Matrices

 
 S = size(dmatnb,1);       % number of spike trains per neuron
 nn= size(dmatnb,2)/S;
 ncmatShift = ones(nn,nn)-eye(nn,nn);  % pre-allocate NxN Neuron Similarity Matrix
 silhHist = zeros(it,crange);  % pre-allocate silhoutte matrix

 
for it = 1:it
    
    randMat = randi(S, [ S nn]);
    
    ix1  = [1:S];
    
      for n1 = 1:nn
        
        refN = dmatnb(randMat(:,n1),ix1)';  % Shuffle rows
        refN = refN(randMat(:,n1),:)';  % Shuffle columns in similar order to rows
       
        ix2  = [1:S];
        
        for n2 = 1:n1
       
            targN = dmatnb(randMat(:,n2),ix2)'; % Shuffle rows
            targN = targN(randMat(:,n2),:)';    % Shuffle columns in similar order to rows
     
            ncmatShift(n1,n2) = corr( refN(:), targN(:) , 'type', corrType);  % Calculate correlation between vectors
            
            if all(refN(:) == 0) & all(targN(:)== 0);       % If BOTH Neurons 'failed' to fire..
           
                ncmatShift(n1,n2) = 1;                      % ...set similarity value to 1 (i.e., the neurons are functionally dissimilar)

            elseif xor([sum(refN(:))],[sum(targN(:))]);      % If either Neuron A or Neuron B exlusively 'failed' to fire...
     
                ncmatShift(n1,n2) = 0;                      % ...set similarity value to 0 (i.e., the neurons are functionally similar)
            end 
            
               ix2 = ix2+S; 
               
        end
  
            ix1 = ix1+S;
   
      end 
    ncmatShift = tril(ncmatShift)+tril(ncmatShift,-1)';         % Copy lower matrix half to upper
    
    % Perform dimensionality Reduction on ncmatShift and calculate
    % silhouette value across the range of tested clusters 
    switch drTechnique
        
        case 'pca'

            pc = pca(ncmatShift);
            NSPACEshift = ncmatShift * pc(:,1:clusterdim);

        case 'mds'

            neuralSpace = cmdscale(1-ncmatShift);   
            NSPACEshift = neuralSpace(:,1:ncmatShift);

        case 'tsne'

            NSPACEshift = runtSNE(ncmatShift,clusterdim,perp); 
    end
    [silhHist(it,:)] = autokmeanscluster(crange,NSPACEshift) ;

end

 shufMU = mean(silhHist,1);     % Calculate the mean silhouette values across all clusters numbers 
 shufCI= prctile(silhHist, [100-silhCI silhCI],1);  % Calculate the upper and low CI (default 0.01)  
  

   
end

