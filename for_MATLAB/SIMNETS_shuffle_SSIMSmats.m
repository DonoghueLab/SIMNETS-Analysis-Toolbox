function [ ncmat, dmatnb, basespiketrains] = SSIMSNeuronRelCorr_mantel(sorted_timestamps, ev, stduration, q)
%%  SPIKE TRAIN SIMILARITY SPACE TOOLBOX% 
% [ncmat] = SSIMSNeuronRelCorr (sorted_timestamps, event, stduration, q)
% helper function for SSIMSETS
%
% @author Carlos Vargas-Irwin
% Copyright (c) Carlos Vargas-Irwin, Brown University. All rightsreserved. Resistance is futile.


%% calculate pairwise distance matrix for each neuron  
[dmatnb, basespiketrains] = getSSIMDMatBetweenTimePoints(sorted_timestamps, ev, stduration, q);
%% Created shuffled ncmat
nn = length(sorted_timestamps); % number of neurons
bl = size(dmatnb,1);
idxTri = tril(ones(bl,bl) , -1); 
ncmat  = zeros(nn,nn);
ncmatS  = zeros(nn,nn);
ix1    = [1:bl];

for ss=1:100

    for n1 = 1:nn

       ix2 = [1:bl];
       for n2 = 1:nn

            if (n1<=n2)

                B = reshape( dmatnb(:,ix2),1,[])';        
                randB = randi(numel(B),[ 1 size(B,1)]);    

                ncmat(n1,n2) = corr( A,B, 'type', 'Pearson');
                ncmatS(n1,n2) = corr( A,B(randB), 'type', 'Pearson');
            else

                ncmat(n1,n2) =   ncmat(n2,n1);
                ncmatS(n1,n2) =   ncmatS(n2,n1);
            end

            ix2 = ix2+bl; 

        end
            ix1 = ix1+bl;
    end 

    corrMat(ss)=corr( ncmat(:),ncmatS(:), 'type', 'Pearson');
end


end

