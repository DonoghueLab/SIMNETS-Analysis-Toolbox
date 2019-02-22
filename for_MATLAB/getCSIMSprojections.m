 function [CSIMSt] = getCSIMSprojections( frmat , metric , dim , perplexity, tSNE_transform, basemat )
% [CSIMS_coordinates, tSNE_transform, basedmat] = getCSIMS(frmat, metric, dim, perplexity)
%
% INPUTS:
% FRmat: matrix of firing rates ( #events x #bins x #neurons)
% metric: used to compare firing rate vectors (see 'help pdist' for options, default = 'euclidean')
% dim: number of desired dimensions for the new space (optional, default = 3)
% perplexity: perplexity value for t-SNE (optional, default = 30)
%
%
% OUTPUT:
% CSIMS_coordinates: dim-dimenstional coordinates of the points associated
% with each firing rate vector
% tSNE_transform: matrix that replicates the transform optimized by tSNE
% 	(multiply new data by this matrix to project onto the same coordinate space)
% basedmat: distance matrix of spike trains. Its size is
% 		[numel(ev), size(frmat,1) * size(frmat,3)]
% 	It consists of the [size(frmat,1), size(frmat,1)] distance matrices for each unit concatenated along the second dimension.

% @author Carlos Vargas-Irwin, Jonas Zimmermann
% Copyright (c) Carlos Vargas-Irwin, Brown University. All rights reserved.

metric = 'euclidean';
dim = 3;
perplexity = 10;
tic;

%dmat = zeros(size(frmat,1), size(frmat,1)*size(frmat,3));


% MATS
dim = size(tSNE_transform, 2);
base_l = size(tSNE_transform, 1);
n_n = size(frmat,3);
n_base_ev = base_l / n_n;
n_ev = size(frmat,1);


dmat = zeros(n_ev, base_l);

% calculate and concatenate distance matrices
ix = 1;
for i_n = 1:size(frmat,3)
    
    
    dmat(:, (1:n_base_ev) + n_base_ev * (i_n - 1) ) = (pdist2(frmat(:,:,i_n),basemat(:,:,i_n), metric));

end

telapsed = toc;
CSIMSt = dmat*tSNE_transform;
 
 




