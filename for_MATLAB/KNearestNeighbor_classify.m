function [knn_class, percorr, vote_matrix, varargout] = KNearestNeighbor_classify(data, true_class, k, varargin)
% [knn_class, percorr, vote_matrix, varargout] = KNearestNeighbor_classify(data, true_class, k, varargin)
% NOTE: data organized in columns (dimensions x number of vectors),
% true_class is a vector of integers representing classes,
% k = # of nearest neighbors to use for classification
% optional 4th parameter: number of permutations for permutation test. Will also have to supply 4th output to
% receive distribution of classification percentages
%
% classifier uses leave-one-out cross validation
%
% @author Carlos Vargas-Irwin, Jonas Zimmermann
% Copyright (c) 2012, 2013, 2014 Carlos Vargas-Irwin, Jonas Zimmermann, Brown University. All rights reserved.

[dim nvec] = size(data);


dmat = pdist(data','euclidean');

%dmat = pdist(data','mahalanobis');
%dmat = pdist(data','seuclidean');
%dmat = pdist(data','chebychev');
%dmat = pdist(data','minkowski');
%dmat = pdist(data','cityblock');


dmat = squareform(dmat);

if nvec <= k
	warning('SSIMS:KNearestNeighbor_classify:k_max', 'In KNearestNeighbor_classify: k = %i, but there are only %i data pints. We will set k to %i.', k, nvec, nvec-1);
	k = nvec - 1;
end

if (nargin > 3) && (nargout > 3)
	n_perm = varargin{1};

	[knn_class, percorr, vote_matrix, perm_class] = knnHelper(dmat, true_class, k, n_perm);

	varargout{1} = perm_class;
else
	[knn_class, percorr, vote_matrix] = knnHelper(dmat, true_class, k);

end

end

function [knn_class, percorr, vote_matrix, perm_class] = knnHelper(dmat, true_class, k, varargin)
	numcorr = 0;
	classes = unique(true_class);
	nvec = size(dmat, 1);

	vote_matrix = zeros(nvec,length(classes));

	perm_class = [];
	n_perm = 0;
	if nargin > 3
		n_perm = varargin{1};
		perm_class = zeros(1, n_perm);
		perm_numcorr = zeros(1, n_perm);

	end

	knn_class = zeros(1, nvec);
	for n = 1:nvec

		% leave one trial out
		nix = [1:n-1 n+1:nvec];
		d = dmat(n,nix);
		%d = d.^2; % squaring shouldn't be necessary - d >= 0 and squaring is strictly monotonous

		[v, ind] = sort(d);

		%randomize ties for first place
		nties = sum(v == v(1));

		[~, ir] = sort(rand(1, nties));
		ti = ind(1:nties);
		ind(1:nties) = ti(ir);

		tc = true_class(nix);
		nnc = tc(ind(1:k));
		[knn_class(n), votes] =  knnFind(classes, nnc);
		if(knn_class(n) == true_class(n))
			numcorr = numcorr + 1;
		end
		vote_matrix(n,:) = votes;
		if n_perm > 0
			for i_p = 1:n_perm
				class_perm = randperm(nvec);
				tc = true_class(class_perm(nix));
				nnc = tc(ind(1:k));
				[p_class, ~] =  knnFind(classes, nnc);
				if(p_class == true_class(class_perm(n)))
					perm_numcorr(i_p) = perm_numcorr(i_p) + 1;
				end
			end
		end
	end
	if n_perm > 0
		perm_class = 100.0 * sort(perm_numcorr) / nvec;
	end
	percorr = 100.0 * ( numcorr / nvec);
end

function [kclass, votes] = knnFind(classes, nnc)
	votes = zeros(numel(classes), 1);
	for c = 1:length(classes)
		votes(c) = sum(nnc == classes(c));
	end

	[~, ci] = max(votes);
	kclass = classes(ci);
end
