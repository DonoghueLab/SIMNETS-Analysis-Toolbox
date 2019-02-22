function [ ncmat, dmatnb, basespiketrains] = SSIMSNeuronRelCorr(sorted_timestamps, ev, stduration, q, varargin)
%% SPIKE TRAIN SIMILARITY SPACE (SSIMS) AND SIMILARITY NETWORKS (SIMNETS) TOOLBOX% 
% OVERVIEW:Steps 1 and 2 of SIMNETS Algorithm
% Step 1: Calculate pairwise distance matrix for each neuron 
% Step 2:Calculate correlation between relational mapping vectors (neurons)
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
% 
% Questions? Contact Carlos_vargas_Irwin@brown.edu or Jacqueline Hynes@Brown.edu
% @author Jacqueline Hynes
% Copyright (c) Jacqueline Hynes, Brown University. All rights reserved.
%% Step 1: Calculate pairwise distance matrix for each neuron 
 
[dmatnb, basespiketrains] = getSSIMDMatBetweenTimePoints(sorted_timestamps, ev, stduration, q);
 

%% Step 2:Calculate correlation between relational mapping vectors (neurons)
numInputs =  numel(varargin); % a and b are required
corrType = {'pearson'}; % placeholder for default values

if ~isempty(varargin{1})
    corrType = varargin{1};
end

   

nn = length(sorted_timestamps);     % number of neurons
S = size(dmatnb,1);     % number of spike trains per neuron
ncmat = zeros(nn,nn);       % pre-allocate NxN Neuron Similarity Matrix
ix1  = [1:S];

% Insert Explanation Here:

for n1 = 1:nn
    
    A = dmatnb(:,ix1);
     
    ix2  = [1:S];
    for n2 = 1:n1
        
        
       B = dmatnb(:,ix2); 

            ncmat(n1,n2) = corr(A(:),B(:) , 'type', corrType);  % Calculate correlation between vectors
  
          % If both neurons are silent, set similarity value to 1 (i.e., the neurons are functionally similar)
          % If either neuron is exclusively silent, set similarity value to 0 (i.e., the neurons are functionally dissimilar)..
            if all(A(:) == 0) & all(B(:)== 0);   % If either Neuron 'failed' to fire, set Corr(A,B)to 1   
           
                ncmat(n1,n2) = 1;      

            elseif xor([sum(A(:))],[sum(B(:))])

                ncmat(n1,n2) = 0;        % If BOTH Neuron A or Neuron B 'failed' to fire, set Corr(A,B)to 1 (i.e., they were both silent)
     
            end 
            

       ix2 = ix2+S; 
        
   end

    ix1 = ix1+S;
               
end 

 ncmat = tril(ncmat)+tril(ncmat,-1)';  % copy lower matrix half to upper  
  
