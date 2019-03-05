function [SSIMS, tSNE_transform, kld] = runtSNE(baseData, outDimensions, perplexity)
% 	RUNTSNE   Perform tSNE dimensionality reduction on high-dimensional data
% 		[SSIMS, TSNE_TRANSFORM, KLD] = RUNTSNE(BASEDATA, OUTDIMENSIONS, PERPLEXITY)
%
% 	@param baseData matrix containing data to be transformed. Rows contain observations, columns contain features
% 	@param outDimensions number of dimensions the dimensionality-reduced data should have. This number should be less than
% 		number of columns of baseData (optional, default = 2)
% 	@param perplexity Sets perplexity parameter for tSNE (optional, default = 30)
%
% 	@return tSNE-transformed data, as well as a tSNE-transformation matrix to transform more data, and also the
% 		Kullback-Leibler-divergence between initial and low-dimensional probability distributions
%
% 	Created by Jonas B Zimmermann on 2016-08-30.
% 	Copyright (c) 2016 Donoghue Lab, Neuroscience Department, Brown University.
% 	All rights reserved.
% 	@author  $Author: Jonas B Zimmermann$
if nargin < 3
    perplexity = 30;
    if nargin < 2
        outDimensions = 3;
    end
end

initial_dims =  size(baseData,1);
%disp('Preprocessing data using PCA...');
if size(baseData, 2) < size(baseData, 1)
    C = baseData'' * baseData;
else
    C = (1 / size(baseData, 1)) * (baseData * baseData');
end
[PCP_transform, lambda] = eig(C);
[lambda, ind] = sort(diag(lambda), 'descend');
PCP_transform = PCP_transform(:,ind(1:initial_dims));
lambda = lambda(1:initial_dims);
if any(lambda <= 0)
    warning('Eigenvalues <= 0! Results may be numerically unstable.\n');
    lambda(lambda<=0) = 1e-8;    % inserted to fix numerical issues
end
if ~(size(baseData, 2) < size(baseData, 1))
    PCP_transform = bsxfun(@times, baseData' * PCP_transform, (1 ./ sqrt(size(baseData, 1) .* lambda))');
end
dmatPCP = baseData * PCP_transform; % no normalization!

[SSIMS, kld] = tsne_PCA(dmatPCP, [], outDimensions, [], perplexity);

%display(['CALCULATING TSNE TRANSFORM ...'])
tSNE_transform = SSIMS' * dmatPCP * pinv(dmatPCP' * dmatPCP);
tSNE_transform = PCP_transform * tSNE_transform';

% added to correct lack of correspondence between SSIMS and transform
% when eigenvalues are corrected (line 55 above) CVI Nov 19 2013
SSIMS = baseData * tSNE_transform;


end %  function

