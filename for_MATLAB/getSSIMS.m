function [SSIMS, tSNE_transform, basespiketrains, basedmat, kld] = getSSIMS(sorted_timestamps, ev, stduration, q, varargin)

% [SSIMS_coordinates] = getSSIMS(sorted_timestamps, ev, stduration, q, dim, perplexity)
%
% INPUTS:
% sorted_timestamps: cell array of neuron spiking stimestamps
% ev: events to align to
% stduration: Length of windows after events to be considered
% q: parameter for spike train metric algorithm (1/q ~ temporal accuracy)
% dim: number of desired dimensions for the new space (optional, default = 3)
% perplexity: perplexity value for t-SNE (optional, default = 30)
%
%
% OUTPUT:
% SSIMS_coordinates: dim-dimenstional coordinates of the points associated with each spike train
% tSNE_transform: matrix that replicates the transform optimized by tSNE
% 	(multiply new data by this matrix to project onto the same coordinate space)
% basespiketrains: spike trains corresponding to events in ev. Cell array,
% 		size(basespiketrains) = [1, numel(sorted_timestamps)]
% 	Each entry in basespiketrains is a cell array of double vectors, and each entry has shape
% 		[1, numel(ev)]
% basedmat: distance matrix of spike trains. Its size is
% 		[numel(ev), numel(ev) * numel(sorted_timestamps)]
% 	It consists of the [numel(ev), numel(ev)] distance matrices for each unit concatenated along the second dimension.

% @author Carlos Vargas-Irwin, Jonas Zimmermann
% Copyright (c) Carlos Vargas-Irwin, Brown University. All rights reserved.
% perplexity = 30;
dim = 3;

warning('SSIMSToolBox:getSSIMS:noCompiledFunction', 'Using plain MATLAB implementation of the SSIMS algorithm. This calculation will be very slow. For speed gains of 1-2 orders of magnitude, consider compilation (see doc/INSTALL.md).');

if ~isempty(varargin)
dim = varargin{1};
end
if nargin > 5
	perplexity = varargin{2};
end

tic;
[basedmat, basespiketrains] = getSSIMDMatBetweenTimePoints(sorted_timestamps, ev, stduration, q);
telapsed = toc;
fprintf('Calculated the base distances for %i units in %4.1f seconds.\n', length(sorted_timestamps), telapsed )

[SSIMS, tSNE_transform, kld] = runtSNE(basedmat, dim, perplexity);

