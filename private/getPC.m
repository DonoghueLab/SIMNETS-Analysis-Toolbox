function [V, ev] = getPC(waves)
% function [V ev] = getPC(waves)
%
% SVD Refresher
% First, make matrix zero-mean
% x = usv'
% u =   (left singular vectors) 
% s(squared) = eigenvalues (diagonal)
% v =  PRINCIPAL COMPONENTS! (right singular vectors) 
%
% note = each column is a sample, each row a variable


zm_wf = zeroMean(waves);


%perform singular value decomposition

[~, S, V] = svd(zm_wf'); 
 
ev = diag(S).^2;