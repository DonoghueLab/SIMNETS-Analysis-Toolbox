 function [CSIMS, tSNE_transform, dmat] = getCSIMS(frmat, varargin)
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
perplexity = 30;


if ~isempty(varargin)
metric = varargin{1};
end
if nargin > 2
	dim = varargin{2};
end
if nargin > 3
	perplexity = varargin{3};
end


tic;
%dmat = zeros(size(frmat,1), size(frmat,1)*size(frmat,3));
% calculate and concatenate distance matrices
ix = 1;
for n = 1:size(frmat,3)
    
   dmat(:,ix:ix+size(frmat,1)-1)= squareform(pdist(frmat(:,:,n),metric));
      %  dmat(:,ix:ix+size(frmat,1)-1)=  (corr(frmat(:,:,n)'));
    
      ix = ix+size(frmat,1);
    
end
telapsed = toc;
% fprintf('Calculated the base distances for %i units in %4.1f seconds.\n', size(frmat,3), telapsed )


% apply t-SNE
[CSIMS, tSNE_transform, kld] = runtSNE(dmat, dim, perplexity);



